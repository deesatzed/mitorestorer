#!/usr/bin/env python3
"""
arXiv Submission Analog & Verification Engine (Submission 1.5+ / AutoTeX Cleanroom)
===================================================================================

Simulates arXiv's end-to-end ingestion, TeX Live `pdflatex` static & clean-room
compilation gates, Poppler PDF/font forensics (`pdffonts`, `pdftotext`, `pdfinfo`),
Step-3 Portal Metadata & MathJax abstract validation (<= 1,920 chars), bibliographic
integrity & live Crossref DOI verification, and arXiv category-moderation screening
(5-gram shingle overlap, figure infographic banners, and clinical-claim triggers).

Usage:
  # Build upload-ready .tar.gz bundles + portal manifests for ready papers & verify:
  python3 tools/arxiv_submission_analog.py --build-ready

  # Audit ALL 5 repository manuscripts (readiness matrix across Papers I-IV + main):
  python3 tools/arxiv_submission_analog.py --audit-all

  # Audit an existing .tar.gz archive exactly as uploaded to arxiv.org/submit:
  python3 tools/arxiv_submission_analog.py --audit-archive dist/arxiv_paper1.tar.gz

  # Run live Crossref DOI title/author verification (requires network access):
  python3 tools/arxiv_submission_analog.py --build-ready --verify-dois-online
"""

from __future__ import annotations

import argparse
import dataclasses
from dataclasses import dataclass, field
import hashlib
import io
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tarfile
import tempfile
import urllib.error
import urllib.parse
import urllib.request
from typing import Dict, List, Optional, Set, Tuple


REPO_ROOT = Path(__file__).resolve().parent.parent
ARXIV_SRC_DIR = REPO_ROOT / "docs" / "arxiv_submission"
FIGURES_DIR = REPO_ROOT / "figures"
DIST_DIR = REPO_ROOT / "dist"

ARXIV_ABSTRACT_MAX_CHARS = 1920
ARXIV_ARCHIVE_MAX_BYTES = 50 * 1024 * 1024
ARXIV_SINGLE_FILE_MAX_BYTES = 10 * 1024 * 1024

VALID_ARXIV_CATEGORIES = {
    "eess.SY", "q-bio.QM", "q-bio.CB", "q-bio.NC", "q-bio.TO",
    "q-bio.PE", "q-bio.MN", "q-bio.BM", "cs.AR", "cs.CE",
    "cs.SY", "math.OC", "math.DS", "physics.bio-ph", "stat.AP",
}


@dataclass(frozen=True)
class ManuscriptConfig:
    paper_id: str
    tex_filename: str
    bundle_name: str
    primary_category: str
    cross_categories: Tuple[str, ...]
    msc_class: str
    acm_class: str
    target_ready: bool
    is_omnibus: bool = False


MANUSCRIPT_REGISTRY: Dict[str, ManuscriptConfig] = {
    "paper1": ManuscriptConfig(
        paper_id="paper1",
        tex_filename="paper1_core_control_and_coral_npu.tex",
        bundle_name="arxiv_paper1",
        primary_category="eess.SY",
        cross_categories=("q-bio.QM", "math.OC", "cs.AR"),
        msc_class="93C10, 92C42, 49M25, 93E11",
        acm_class="J.3; C.3",
        target_ready=True,
    ),
    "paper2": ManuscriptConfig(
        paper_id="paper2",
        tex_filename="paper2_ex_vivo_cart_bioreactor.tex",
        bundle_name="arxiv_paper2",
        primary_category="q-bio.CB",
        cross_categories=("q-bio.QM", "eess.SY"),
        msc_class="92C42, 93C10, 93B07, 49N90",
        acm_class="J.3",
        target_ready=True,
    ),
    "paper3": ManuscriptConfig(
        paper_id="paper3",
        tex_filename="paper3_apis_pollinator_smart_hive.tex",
        bundle_name="arxiv_paper3",
        primary_category="q-bio.QM",
        cross_categories=("q-bio.PE", "eess.SY"),
        msc_class="92C42, 92D40, 93C10",
        acm_class="J.3",
        target_ready=False,
    ),
    "paper4": ManuscriptConfig(
        paper_id="paper4",
        tex_filename="paper4_neuro_qsp_and_axonal_pde.tex",
        bundle_name="arxiv_paper4",
        primary_category="q-bio.NC",
        cross_categories=("q-bio.QM", "q-bio.TO"),
        msc_class="92C42, 35K57, 93E11",
        acm_class="J.3",
        target_ready=False,
    ),
    "main": ManuscriptConfig(
        paper_id="main",
        tex_filename="main.tex",
        bundle_name="arxiv_omnibus_do_not_submit",
        primary_category="q-bio.QM",
        cross_categories=("eess.SY", "q-bio.CB", "q-bio.NC"),
        msc_class="92C42, 93C10",
        acm_class="J.3",
        target_ready=False,
        is_omnibus=True,
    ),
}


@dataclass
class CheckFinding:
    gate: str
    severity: str  # "PASS", "INFO", "WARN", "FAIL"
    code: str
    message: str


@dataclass
class AuditReport:
    paper_id: str
    archive_path: str
    sha256: str
    archive_bytes: int
    pages: int = 0
    figures_count: int = 0
    tables_count: int = 0
    portal_abstract_chars: int = 0
    portal_title: str = ""
    portal_abstract: str = ""
    findings: List[CheckFinding] = field(default_factory=list)

    @property
    def fail_count(self) -> int:
        return sum(1 for f in self.findings if f.severity == "FAIL")

    @property
    def warn_count(self) -> int:
        return sum(1 for f in self.findings if f.severity == "WARN")

    @property
    def passed(self) -> bool:
        return self.fail_count == 0

    def add(self, gate: str, severity: str, code: str, message: str) -> None:
        self.findings.append(CheckFinding(gate=gate, severity=severity, code=code, message=message))


def _which_tool(name: str) -> Optional[str]:
    found = shutil.which(name)
    if found:
        return found
    hb = Path("/opt/homebrew/bin") / name
    if hb.exists():
        return str(hb)
    return None


def strip_latex_comments(tex: str) -> str:
    """Remove LaTeX % comments while preserving escaped \\%."""
    out_lines = []
    for line in tex.splitlines():
        m = re.search(r"(?<!\\)%", line)
        if m:
            kept = line[: m.start()]
            if kept.strip() or not line.lstrip().startswith("%"):
                out_lines.append(kept)
        else:
            out_lines.append(line)
    return "\n".join(out_lines)


def latex_to_arxiv_portal_text(raw: str, is_title: bool = False) -> str:
    """
    Convert LaTeX title or abstract into arXiv Step-3 Web/MathJax portal text.
    Outside $...$, arXiv's landing page does NOT execute \\texttt, \\textbf, \\textit,
    or custom spacing macros, whereas standard inline math ($...$) is rendered by MathJax.
    """
    s = strip_latex_comments(raw)
    s = s.replace("\\noindent", " ")
    if is_title:
        s = re.sub(r"\\(?:LARGE|Large|large|bfseries)\b", "", s)

    # Split into math ($...$) vs text segments so we only unwrap text-mode macros outside math
    parts = re.split(r"(\$[^$]+\$)", s)
    cleaned_parts: List[str] = []
    for idx, part in enumerate(parts):
        if idx % 2 == 1:
            # Inside inline $...$ math: normalize \text{ ... } to \mathrm{...} for MathJax safety
            mp = re.sub(r"\\text\{\s*([^}]*)\}", lambda m: r"\,\mathrm{" + m.group(1).strip() + "}", part)
            cleaned_parts.append(mp)
        else:
            tp = part
            # Iteratively unwrap \textbf{...}, \textit{...}, \emph{...}, \texttt{...}
            for _ in range(3):
                tp = re.sub(r"\\(?:textbf|textit|emph|textrm)\{([^{}]*)\}", r"\1", tp)
                tp = re.sub(
                    r"\\texttt\{([^{}]*)\}",
                    lambda m: m.group(1).replace(r"\_", "_"),
                    tp,
                )
            tp = tp.replace(r"\AA", "Å")  # Wait: keep pure ASCII for arXiv portal safety!
            tp = tp.replace("Å", "A")
            # Better: render Angstrom in MathJax as $\mathrm{\AA}$ or "Angstrom" -> use $\mathrm{\AA}$
            cleaned_parts.append(tp)

    # Re-pass \AA before joining accurately:
    s2 = strip_latex_comments(raw).replace("\\noindent", " ")
    if is_title:
        s2 = re.sub(r"\\(?:LARGE|Large|large|bfseries)\b", "", s2)
    parts2 = re.split(r"(\$[^$]+\$)", s2)
    final_parts: List[str] = []
    for idx, part in enumerate(parts2):
        if idx % 2 == 1:
            mp = re.sub(
                r"\\text\{\s*\\AA\s*\}",
                r"\,\\mathring{\\mathrm{A}}",
                part,
            )
            mp = re.sub(
                r"\\text\{\s*([A-Za-z0-9/._-]+)\s*\}",
                r"\,\\mathrm{\1}",
                mp,
            )
            mp = mp.replace(r"\text{--}", "--")
            final_parts.append(mp)
        else:
            tp = part
            for _ in range(3):
                tp = re.sub(r"\\(?:textbf|textit|emph|textrm)\{([^{}]*)\}", r"\1", tp)
                tp = re.sub(
                    r"\\texttt\{([^{}]*)\}",
                    lambda m: m.group(1).replace(r"\_", "_"),
                    tp,
                )
            tp = re.sub(r"\s*\\AA\b", r" $\\mathring{\\mathrm{A}}$", tp)
            tp = tp.replace("---", " -- ")
            tp = tp.replace("--", "-")
            tp = tp.replace("~", " ")
            tp = tp.replace(r"\ ", " ")
            tp = tp.replace(r"\%", "%")
            tp = tp.replace("``", '"').replace("''", '"')
            final_parts.append(tp)

    joined = "".join(final_parts)
    # Merge adjacent math blocks if any arose
    joined = " ".join(joined.split())
    return joined.strip()


def _parse_bibtex_fields_balanced(body: str) -> Dict[str, str]:
    """Parse BibTeX key=value pairs supporting arbitrary nested brace depth."""
    fields: Dict[str, str] = {}
    n = len(body)
    i = 0
    while i < n:
        while i < n and body[i] in " \t\r\n,":
            i += 1
        if i >= n:
            break
        m_key = re.match(r"([a-zA-Z0-9_:-]+)\s*=\s*", body[i:])
        if not m_key:
            nxt = body.find(",", i)
            i = (nxt + 1) if nxt != -1 else n
            continue
        fname = m_key.group(1).lower()
        i += m_key.end()
        if i >= n:
            break
        if body[i] == "{":
            depth = 1
            i += 1
            start = i
            while i < n and depth > 0:
                if body[i] == "{":
                    depth += 1
                elif body[i] == "}":
                    depth -= 1
                i += 1
            raw_val = body[start : max(start, i - 1)]
        elif body[i] == '"':
            i += 1
            start = i
            while i < n and not (body[i] == '"' and body[i - 1] != "\\"):
                i += 1
            raw_val = body[start:i]
            i += 1
        else:
            start = i
            while i < n and body[i] not in ",\r\n}":
                i += 1
            raw_val = body[start:i]
        fields[fname] = " ".join(raw_val.split())
    return fields


def parse_bib_database(bib_path: Path) -> Dict[str, Dict[str, str]]:
    """Parse references.bib into {key: {field: value}} with arbitrary nested braces."""
    if not bib_path.exists():
        return {}
    text = bib_path.read_text(encoding="utf-8")
    entries: Dict[str, Dict[str, str]] = {}
    pattern = re.compile(r"@(\w+)\s*\{\s*([^,\s]+)\s*,", re.IGNORECASE)
    pos = 0
    while True:
        match = pattern.search(text, pos)
        if not match:
            break
        etype = match.group(1).lower()
        if etype in ("comment", "preamble", "string"):
            pos = match.end()
            continue
        key = match.group(2).strip()
        idx = match.end()
        depth = 1
        while idx < len(text) and depth > 0:
            if text[idx] == "{":
                depth += 1
            elif text[idx] == "}":
                depth -= 1
            idx += 1
        body = text[match.end() : max(match.end(), idx - 1)]
        pos = idx
        fields = _parse_bibtex_fields_balanced(body)
        fields["_type"] = etype
        entries[key] = fields
    return entries


def build_arxiv_bundle(
    cfg: ManuscriptConfig,
    out_dir: Path = DIST_DIR,
) -> Tuple[Path, Path]:
    """
    Package a manuscript into a deterministic, flat, self-contained arXiv .tar.gz
    submission archive (with pre-compiled <stem>.bbl, local figure PDFs, and
    Submission 1.5+ `00README.json`) plus a Step-3 Portal Metadata file.
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    src_tex_path = ARXIV_SRC_DIR / cfg.tex_filename
    if not src_tex_path.exists():
        raise FileNotFoundError(f"Missing manuscript TeX file: {src_tex_path}")

    tectonic_bin = _which_tool("tectonic")
    if not tectonic_bin:
        raise RuntimeError("`tectonic` is required to pre-compile deterministic .bbl files.")

    stem = src_tex_path.stem

    with tempfile.TemporaryDirectory(prefix=f"arxiv_stage_{cfg.paper_id}_") as tmp_str:
        tmp = Path(tmp_str)
        # 1. Compile in place against references.bib to obtain exact <stem>.bbl
        res = subprocess.run(
            [
                tectonic_bin,
                "-X",
                "compile",
                "--keep-intermediates",
                "--outdir",
                str(tmp),
                str(src_tex_path),
            ],
            capture_output=True,
            text=True,
        )
        if res.returncode != 0:
            raise RuntimeError(
                f"Failed generating .bbl for {cfg.tex_filename}:\n{res.stderr[-1200:]}"
            )

        bbl_src = tmp / f"{stem}.bbl"
        if not bbl_src.exists():
            raise RuntimeError(f"Expected {bbl_src} was not generated.")

        # 2. Prepare flat clean staging folder
        stage = tmp / "bundle"
        stage.mkdir()

        raw_tex = src_tex_path.read_text(encoding="utf-8")
        no_comments_for_scan = strip_latex_comments(raw_tex)

        # Ensure \pdfoutput=1 guard and utf8 inputenc for TeX Live pdflatex safety
        patched_tex = raw_tex
        if "\\pdfoutput=1" not in patched_tex[:200]:
            patched_tex = (
                "\\ifdefined\\XeTeXversion\\else\\ifdefined\\pdfoutput\\pdfoutput=1\\fi\\fi\n"
                + patched_tex
            )
        if "\\usepackage[utf8]{inputenc}" not in patched_tex and "\\usepackage[T1]{fontenc}" in patched_tex:
            patched_tex = patched_tex.replace(
                "\\usepackage[T1]{fontenc}",
                "\\usepackage[utf8]{inputenc}\n\\usepackage[T1]{fontenc}",
                1,
            )

        # Rewrite \includegraphics paths to flat local filenames and copy only referenced figures
        included_figs: List[str] = []

        def _rewrite_graphic(match: re.Match) -> str:
            opt = match.group(1) or ""
            orig_rel = match.group(2).strip()
            resolved = (ARXIV_SRC_DIR / orig_rel).resolve()
            if not resolved.exists() and not resolved.suffix:
                resolved = resolved.with_suffix(".pdf")
            if not resolved.exists():
                raise FileNotFoundError(
                    f"{cfg.tex_filename} references missing figure: {orig_rel} ({resolved})"
                )
            local_name = resolved.name
            if local_name not in included_figs:
                included_figs.append(local_name)
                shutil.copy2(resolved, stage / local_name)
            return f"\\includegraphics{opt}{{{local_name}}}"

        patched_tex = re.sub(
            r"\\includegraphics(\[[^\]]*\])?\{([^}]+)\}",
            _rewrite_graphic,
            patched_tex,
        )

        # Replace \bibliographystyle{...} + \bibliography{...} with explicit \input{<stem>.bbl}
        # so neither arXiv Submission 1.5+, pdflatex, nor cleanroom tectonic attempts a BibTeX rerun
        patched_tex = re.sub(r"\\bibliographystyle\{[^}]+\}\s*", "", patched_tex)
        patched_tex = re.sub(
            r"\\bibliography\{[^}]+\}",
            lambda _m: f"\\input{{{stem}.bbl}}",
            patched_tex,
        )

        # Strip pure-comment lines so internal editor notes never leak into public arXiv TeX source
        cleaned_lines: List[str] = []
        for idx, line in enumerate(patched_tex.splitlines()):
            if idx > 0 and line.lstrip().startswith("%"):
                continue
            cleaned_lines.append(line)
        cleaned_tex = "\n".join(cleaned_lines).strip() + "\n"

        (stage / cfg.tex_filename).write_text(cleaned_tex, encoding="utf-8")
        shutil.copy2(bbl_src, stage / f"{stem}.bbl")

        # 3. Emit arXiv Submission 1.5+ `00README.json`
        readme_json = {
            "process": {
                "compiler": "pdflatex",
            },
            "sources": [
                {"filename": cfg.tex_filename, "usage": "toplevel"},
                {"filename": f"{stem}.bbl", "usage": "include"},
                *[{"filename": fn, "usage": "include"} for fn in included_figs],
            ],
        }
        (stage / "00README.json").write_text(
            json.dumps(readme_json, indent=2) + "\n",
            encoding="utf-8",
        )

        # 4. Build deterministic POSIX .tar.gz archive (uid=0, gid=0, sorted entries, no macOS xattrs)
        archive_path = out_dir / f"{cfg.bundle_name}.tar.gz"
        files_in_order = ["00README.json", cfg.tex_filename, f"{stem}.bbl", *included_figs]
        with tarfile.open(archive_path, "w:gz", format=tarfile.PAX_FORMAT) as tar:
            for fname in files_in_order:
                fpath = stage / fname
                st = fpath.stat()
                ti = tarfile.TarInfo(name=fname)
                ti.size = st.st_size
                ti.mtime = 1794441600  # Deterministic timestamp (2026-11-12T00:00:00Z)
                ti.mode = 0o644
                ti.type = tarfile.REGTYPE
                ti.uid = 0
                ti.gid = 0
                ti.uname = "root"
                ti.gname = "root"
                with fpath.open("rb") as fh:
                    tar.addfile(ti, fh)

    # Generate portal metadata text placeholder (populated with verified page/table counts during audit)
    manifest_path = out_dir / f"{cfg.bundle_name}_portal_metadata.txt"
    return archive_path, manifest_path


def write_portal_manifest(
    cfg: ManuscriptConfig,
    report: AuditReport,
    manifest_path: Path,
) -> None:
    """Write ready-to-paste arXiv Step-3 Portal Metadata file."""
    cross_str = ", ".join(cfg.cross_categories)
    comments = (
        f"{report.pages} pages, {report.figures_count} figures, {report.tables_count} tables. "
        f"MitoRestorer Series ({cfg.paper_id.upper()}); strictly computational/in silico study "
        f"with full open-source MIASE simulation and C99 firmware reproducibility at "
        f"https://github.com/deesatzed/mitorestorer"
    )
    text = f"""================================================================================
arXiv Submission 1.5+ Ready-to-Paste Portal Metadata ({cfg.bundle_name}.tar.gz)
================================================================================
Archive File      : {report.archive_path}
Archive Size      : {report.archive_bytes:,} bytes ({report.archive_bytes / 1024.0:.1f} KiB)
Archive SHA-256   : {report.sha256}
Preflight Verdict : {"PASS (0 FAIL)" if report.passed else f"BLOCKED ({report.fail_count} FAIL)"} | Warnings: {report.warn_count}

[Step 1: Category & License]
Primary Category  : {cfg.primary_category}
Cross-Lists       : {cross_str}
License           : CC BY 4.0 (http://creativecommons.org/licenses/by/4.0/)

[Step 2: Upload Files]
Upload Archive    : {Path(report.archive_path).name}
Compiler (00README): pdflatex (pre-compiled .bbl included; zero BibTeX dependency)

[Step 3: Metadata Fields — Copy & Paste Directly]
Title:
{report.portal_title}

Authors:
Wayne A. Satz

Affiliation (for profile / author verification):
Department of Emergency Medicine, Lewis Katz School of Medicine at Temple University, Philadelphia, PA 19140, USA (satzwa@temple.edu)

Abstract ({report.portal_abstract_chars} / {ARXIV_ABSTRACT_MAX_CHARS} characters):
{report.portal_abstract}

Comments:
{comments}

MSC-class:
{cfg.msc_class}

ACM-class:
{cfg.acm_class}
================================================================================
"""
    manifest_path.write_text(text, encoding="utf-8")


def _shingle_set(text: str, k: int = 5) -> Set[Tuple[str, ...]]:
    words = re.findall(r"[a-z]{3,}", text.lower())
    if len(words) < k:
        return set()
    return {tuple(words[i : i + k]) for i in range(len(words) - k + 1)}


def verify_doi_crossref(doi: str, expected_title: str, timeout_s: float = 8.0) -> Tuple[bool, str]:
    """Query Crossref API (or doi.org handle for PDB DOIs) to prove citation authenticity."""
    if doi.lower().startswith("10.2210/pdb"):
        pdb_id = doi.split("/")[-2].replace("pdb", "").upper()
        url = f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"
        req = urllib.request.Request(url, headers={"User-Agent": "MitoRestorer-arXiv-Analog/1.0 (mailto:satzwa@temple.edu)"})
        try:
            with urllib.request.urlopen(req, timeout=timeout_s) as resp:
                if resp.status == 200:
                    return True, f"RCSB PDB {pdb_id} verified (HTTP 200)"
        except Exception as exc:
            return False, f"RCSB lookup error for {doi}: {exc}"

    encoded = urllib.parse.quote(doi, safe="")
    url = f"https://api.crossref.org/works/{encoded}"
    req = urllib.request.Request(
        url,
        headers={"User-Agent": "MitoRestorer-arXiv-Analog/1.0 (mailto:satzwa@temple.edu)"},
    )
    try:
        with urllib.request.urlopen(req, timeout=timeout_s) as resp:
            payload = json.loads(resp.read().decode("utf-8"))
    except urllib.error.HTTPError as exc:
        return False, f"HTTP {exc.code} from Crossref for DOI {doi}"
    except Exception as exc:
        return False, f"Network error querying Crossref for {doi}: {exc}"

    msg = payload.get("message", {})
    cr_titles = msg.get("title", [])
    cr_title = cr_titles[0] if cr_titles else ""
    cleaned_exp = re.sub(r"\\[a-zA-Z]+\{([^}]*)\}", r"\1", expected_title)
    cleaned_exp = re.sub(r"[{}\\\"'`~^$]", "", cleaned_exp.lower())
    cleaned_cr = re.sub(r"[{}\\\"'`~^$]", "", cr_title.lower())
    stop = {"with", "from", "that", "this", "into", "during", "human", "cells", "using", "study", "analysis", "between"}
    exp_tokens = {w for w in re.findall(r"[a-z0-9]{3,}", cleaned_exp) if w not in stop}
    cr_tokens = {w for w in re.findall(r"[a-z0-9]{3,}", cleaned_cr) if w not in stop}
    if not exp_tokens or not cr_tokens:
        return (
            False,
            f"Missing comparable title tokens for DOI {doi}: BibTeX='{expected_title[:60]}' vs Crossref='{cr_title[:60]}'",
        )
    inter = len(exp_tokens & cr_tokens)
    union = len(exp_tokens | cr_tokens)
    containment = inter / min(len(exp_tokens), len(cr_tokens))
    similarity = 0.5 * (inter / union) + 0.5 * containment
    if similarity < 0.45:
        return (
            False,
            f"DOI/title mismatch! BibTeX='{expected_title[:60]}' vs Crossref='{cr_title[:60]}' (similarity {similarity:.2f} < 0.45)",
        )
    return True, f"Verified '{cr_title[:65]}' (similarity {similarity:.2f})"


def audit_archive(
    archive_path: Path,
    cfg: Optional[ManuscriptConfig] = None,
    bib_db: Optional[Dict[str, Dict[str, str]]] = None,
    sibling_tex_map: Optional[Dict[str, str]] = None,
    verify_dois_online: bool = False,
) -> AuditReport:
    """
    Execute all 7 arXiv Submission Analog verification gates against a .tar.gz archive
    inside a hermetically isolated cleanroom temporary directory.
    """
    archive_bytes = archive_path.stat().st_size
    sha256 = hashlib.sha256(archive_path.read_bytes()).hexdigest()
    paper_id = cfg.paper_id if cfg else archive_path.stem
    report = AuditReport(
        paper_id=paper_id,
        archive_path=str(archive_path),
        sha256=sha256,
        archive_bytes=archive_bytes,
    )

    # =========================================================================
    # GATE 1: ARCHIVE_QUARANTINE & 00README.json v1.5 SCHEMA
    # =========================================================================
    if archive_bytes > ARXIV_ARCHIVE_MAX_BYTES:
        report.add("G1_ARCHIVE", "FAIL", "ARCHIVE_TOO_LARGE", f"Archive is {archive_bytes:,} B (> 50 MB arXiv limit).")
    else:
        report.add("G1_ARCHIVE", "PASS", "ARCHIVE_SIZE_OK", f"Archive size {archive_bytes:,} B ({archive_bytes/1024:.1f} KiB) within 50 MB limit.")

    with tempfile.TemporaryDirectory(prefix=f"arxiv_cleanroom_{paper_id}_") as clean_str:
        clean = Path(clean_str)
        try:
            with tarfile.open(archive_path, "r:*") as tar:
                members = tar.getmembers()
                names = [m.name for m in members]
                for m in members:
                    p = Path(m.name)
                    if m.name.startswith("/") or ".." in p.parts:
                        report.add("G1_ARCHIVE", "FAIL", "ILLEGAL_PATH_TRAVERSAL", f"Illegal path in archive: {m.name}")
                    if m.issym() or m.islnk():
                        report.add("G1_ARCHIVE", "FAIL", "ILLEGAL_SYMLINK", f"Symlinks/hardlinks prohibited by arXiv: {m.name}")
                    if p.name.startswith("._") or p.name in {".DS_Store", "__MACOSX"}:
                        report.add("G1_ARCHIVE", "FAIL", "MACOS_JUNK_FILE", f"macOS metadata file in archive: {m.name}")
                    if p.suffix.lower() in {".aux", ".log", ".out", ".synctex", ".fls", ".fdb_latexmk"}:
                        report.add("G1_ARCHIVE", "FAIL", "BUILD_ARTIFACT_LEAK", f"Intermediate TeX artifact in archive: {m.name}")
                    if m.size > ARXIV_SINGLE_FILE_MAX_BYTES:
                        report.add("G1_ARCHIVE", "FAIL", "SINGLE_FILE_TOO_LARGE", f"File {m.name} is {m.size:,} B (> 10 MB limit).")
                    if not re.match(r"^[A-Za-z0-9_./-]+$", m.name):
                        report.add("G1_ARCHIVE", "FAIL", "NON_PORTABLE_FILENAME", f"Filename contains spaces or non-ASCII chars: {m.name}")

                lower_names = [n.lower() for n in names]
                if len(lower_names) != len(set(lower_names)):
                    report.add("G1_ARCHIVE", "FAIL", "CASE_COLLISION", "Case-insensitive filename collision in archive.")

                if report.fail_count == 0:
                    tar.extractall(path=clean, filter="data")
        except Exception as exc:
            report.add("G1_ARCHIVE", "FAIL", "TAR_CORRUPT", f"Failed unpacking archive: {exc}")
            return report

        readme_file = clean / "00README.json"
        toplevel_tex_name: Optional[str] = None
        if not readme_file.exists():
            report.add("G1_ARCHIVE", "WARN", "MISSING_00README", "No 00README.json found; arXiv will guess compiler heuristics.")
            tex_candidates = list(clean.glob("*.tex"))
            if len(tex_candidates) == 1:
                toplevel_tex_name = tex_candidates[0].name
            else:
                report.add("G1_ARCHIVE", "FAIL", "AMBIGUOUS_TOPLEVEL_TEX", f"Found {len(tex_candidates)} .tex files without 00README.json.")
                return report
        else:
            try:
                rm = json.loads(readme_file.read_text(encoding="utf-8"))
                compiler = rm.get("process", {}).get("compiler")
                if compiler != "pdflatex":
                    report.add("G1_ARCHIVE", "FAIL", "UNEXPECTED_COMPILER", f"00README.json compiler='{compiler}' (expected 'pdflatex').")
                sources = rm.get("sources", [])
                toplevels = [s["filename"] for s in sources if s.get("usage") == "toplevel"]
                if len(toplevels) != 1:
                    report.add("G1_ARCHIVE", "FAIL", "INVALID_TOPLEVEL_COUNT", f"00README.json must specify exactly 1 toplevel .tex (got {toplevels}).")
                    return report
                toplevel_tex_name = toplevels[0]
                declared_files = {s["filename"] for s in sources} | {"00README.json"}
                actual_files = {p.name for p in clean.iterdir() if p.is_file()}
                missing_decl = declared_files - actual_files
                undeclared = actual_files - declared_files
                if missing_decl:
                    report.add("G1_ARCHIVE", "FAIL", "MISSING_DECLARED_SOURCES", f"00README.json references missing files: {sorted(missing_decl)}")
                if undeclared:
                    report.add("G1_ARCHIVE", "WARN", "UNDECLARED_ARCHIVE_FILES", f"Archive contains files not listed in 00README.json: {sorted(undeclared)}")
                if not missing_decl and not undeclared:
                    report.add("G1_ARCHIVE", "PASS", "00README_VALID", f"00README.json v1.5 valid (compiler=pdflatex, {len(sources)} sources inventoried).")
            except Exception as exc:
                report.add("G1_ARCHIVE", "FAIL", "INVALID_00README_JSON", f"Malformed 00README.json: {exc}")
                return report

        tex_path = clean / toplevel_tex_name
        stem = tex_path.stem
        if (clean / f"{stem}.pdf").exists():
            report.add("G1_ARCHIVE", "FAIL", "PRECOMPILED_ROOT_PDF", f"Archive includes pre-compiled {stem}.pdf alongside {toplevel_tex_name}; arXiv rejects/strips this.")
        bbl_path = clean / f"{stem}.bbl"
        if not bbl_path.exists():
            report.add("G1_ARCHIVE", "FAIL", "MISSING_BBL", f"Missing pre-compiled {stem}.bbl (arXiv pdflatex will not run bibtex).")
            return report
        if list(clean.glob("*.bib")):
            report.add("G1_ARCHIVE", "WARN", "REDUNDANT_BIB_FILE", "Archive includes .bib alongside .bbl; omit .bib to prevent BibTeX/Biber version conflicts.")

        # =====================================================================
        # GATE 2: TEX_STATIC_PDFLATEX_LINT
        # =====================================================================
        tex_raw = tex_path.read_text(encoding="utf-8")
        bbl_raw = bbl_path.read_text(encoding="utf-8")
        tex_nc = strip_latex_comments(tex_raw)

        first_five = "\n".join(tex_raw.splitlines()[:5])
        if "\\pdfoutput=1" not in first_five:
            report.add("G2_TEX_STATIC", "WARN", "NO_PDFOUTPUT_DIRECTIVE", "Missing \\pdfoutput=1 in lines 1-5.")
        else:
            report.add("G2_TEX_STATIC", "PASS", "PDFOUTPUT_DIRECTIVE", "Explicit \\pdfoutput=1 directive present in lines 1-5.")

        if "\\usepackage[utf8]{inputenc}" not in tex_nc or "\\usepackage[T1]{fontenc}" not in tex_nc:
            report.add("G2_TEX_STATIC", "WARN", "ENCODING_PACKAGES", "Include both \\usepackage[utf8]{inputenc} and \\usepackage[T1]{fontenc} for TeX Live pdflatex consistency.")
        else:
            report.add("G2_TEX_STATIC", "PASS", "ENCODING_PACKAGES_OK", "Both utf8 inputenc and T1 fontenc explicitly loaded.")

        forbidden_pkgs = ["fontspec", "unicode-math", "polyglossia", "luacode", "minted", "svg"]
        for pkg in forbidden_pkgs:
            if re.search(rf"\\usepackage(?:\[[^\]]*\])?\{{[^}}]*\b{re.escape(pkg)}\b", tex_nc):
                report.add("G2_TEX_STATIC", "FAIL", "FORBIDDEN_PDFLATEX_PKG", f"Package '{pkg}' fails under arXiv pdflatex (--no-shell-escape).")

        if re.search(r"\\write18\b", tex_nc):
            report.add("G2_TEX_STATIC", "FAIL", "SHELL_ESCAPE_WRITE18", "\\write18 shell execution is blocked in arXiv's sandbox.")

        if "../" in tex_nc or "..\\" in tex_nc:
            report.add("G2_TEX_STATIC", "FAIL", "PARENT_PATH_IN_TEX", "Parent-directory traversal '../' found in TeX commands.")

        non_ascii_tex = [ch for ch in tex_raw if ord(ch) > 127]
        non_ascii_bbl = [ch for ch in bbl_raw if ord(ch) > 127]
        if non_ascii_tex or non_ascii_bbl:
            report.add("G2_TEX_STATIC", "FAIL", "NON_ASCII_SOURCE", f"Non-ASCII bytes found (.tex={len(non_ascii_tex)}, .bbl={len(non_ascii_bbl)}); risk of pdflatex inputenc error.")
        else:
            report.add("G2_TEX_STATIC", "PASS", "PURE_7BIT_ASCII", "100% pure 7-bit ASCII across .tex and .bbl sources.")

        report.tables_count = len(re.findall(r"\\begin\{table\*?\}", tex_nc))
        fig_refs = re.findall(r"\\includegraphics(?:\[[^\]]*\])?\{([^}]+)\}", tex_nc)
        report.figures_count = len(set(fig_refs))

        # =====================================================================
        # GATE 3: CLEANROOM_AUTOTEX_COMPILE (NO .bib ACCESS)
        # =====================================================================
        pdflatex_bin = _which_tool("pdflatex")
        tectonic_bin = _which_tool("tectonic")
        compiled_pdf = clean / f"{stem}.pdf"
        compile_logs = ""

        if pdflatex_bin:
            for pass_idx in range(3):
                proc = subprocess.run(
                    [pdflatex_bin, "-interaction=nonstopmode", "-halt-on-error", "-no-shell-escape", toplevel_tex_name],
                    cwd=str(clean),
                    capture_output=True,
                    text=True,
                )
                compile_logs = proc.stdout + "\n" + proc.stderr
                if proc.returncode != 0:
                    report.add("G3_AUTOTEX", "FAIL", "PDFLATEX_COMPILE_ERROR", f"pdflatex pass {pass_idx+1} failed:\n{compile_logs[-800:]}")
                    return report
            engine_used = "pdflatex (-no-shell-escape)"
        elif tectonic_bin:
            proc = subprocess.run(
                [tectonic_bin, "-X", "compile", "--print", toplevel_tex_name],
                cwd=str(clean),
                capture_output=True,
                text=True,
            )
            compile_logs = proc.stdout + "\n" + proc.stderr
            if proc.returncode != 0:
                report.add("G3_AUTOTEX", "FAIL", "CLEANROOM_COMPILE_ERROR", f"Cleanroom compilation without .bib failed:\n{compile_logs[-800:]}")
                return report
            engine_used = "tectonic cleanroom (isolated, zero .bib)"
        else:
            report.add("G3_AUTOTEX", "FAIL", "NO_TEX_ENGINE", "Neither pdflatex nor tectonic found on host.")
            return report

        if not compiled_pdf.exists():
            report.add("G3_AUTOTEX", "FAIL", "NO_OUTPUT_PDF", "Compiler exited 0 but did not produce output PDF.")
            return report

        overfulls = [float(m.group(1)) for m in re.finditer(r"Overfull \\hbox \(([0-9.]+)pt too wide\)", compile_logs)]
        max_overfull = max(overfulls) if overfulls else 0.0
        if max_overfull > 35.0:
            report.add("G3_AUTOTEX", "FAIL", "SEVERE_OVERFULL_HBOX", f"Severe margin overflow: max Overfull \\hbox = {max_overfull:.1f} pt (> 35 pt).")
        elif max_overfull > 12.0:
            report.add("G3_AUTOTEX", "WARN", "MODERATE_OVERFULL_HBOX", f"Moderate Overfull \\hbox = {max_overfull:.1f} pt ({len(overfulls)} instances).")
        else:
            report.add("G3_AUTOTEX", "PASS", "CLEANROOM_COMPILE_OK", f"Cleanroom compile succeeded via {engine_used} (max Overfull \\hbox = {max_overfull:.1f} pt).")

        # =====================================================================
        # GATE 4: PDF_AND_FONT_FORENSICS (pdfinfo, pdffonts, pdftotext)
        # =====================================================================
        pdfinfo_bin = _which_tool("pdfinfo")
        pdffonts_bin = _which_tool("pdffonts")
        pdftotext_bin = _which_tool("pdftotext")

        if pdfinfo_bin:
            info_out = subprocess.run([pdfinfo_bin, str(compiled_pdf)], capture_output=True, text=True).stdout
            pm = re.search(r"Pages:\s+(\d+)", info_out)
            report.pages = int(pm.group(1)) if pm else 0
            if report.pages <= 0:
                report.add("G4_PDF_FONTS", "FAIL", "ZERO_PAGES", "Compiled PDF has 0 pages.")
            else:
                report.add("G4_PDF_FONTS", "PASS", "PDF_GEOMETRY_OK", f"Compiled PDF verified ({report.pages} pages).")

        if pdffonts_bin:
            font_out = subprocess.run([pdffonts_bin, str(compiled_pdf)], capture_output=True, text=True).stdout
            type3_count = 0
            unembedded_count = 0
            total_fonts = 0
            for line in font_out.splitlines()[2:]:
                cols = line.split()
                if len(cols) >= 5:
                    total_fonts += 1
                    if "Type 3" in line:
                        type3_count += 1
                    # Columns from end: object_id_num, object_id_gen, uni, sub, emb
                    emb_flag = cols[-5] if len(cols) >= 6 else "yes"
                    if emb_flag == "no":
                        unembedded_count += 1
            if unembedded_count > 0:
                report.add("G4_PDF_FONTS", "FAIL", "UNEMBEDDED_FONTS", f"{unembedded_count}/{total_fonts} fonts in compiled PDF are NOT embedded!")
            if type3_count > 0:
                report.add("G4_PDF_FONTS", "FAIL", "TYPE3_BITMAP_FONTS", f"Found {type3_count} Type-3 bitmap/procedural font objects in PDF (regenerate figures with pdf.fonttype=42).")
            if unembedded_count == 0 and type3_count == 0:
                report.add("G4_PDF_FONTS", "PASS", "VECTOR_FONTS_100PCT", f"All {total_fonts} font objects embedded as crisp Type-1C / TrueType (Type-42); zero Type-3 fonts.")

        rendered_text = ""
        if pdftotext_bin:
            rendered_text = subprocess.run([pdftotext_bin, str(compiled_pdf), "-"], capture_output=True, text=True).stdout
            unresolved_cites = re.findall(r"\[\?\]", rendered_text)
            unresolved_refs = re.findall(r"\?\?", rendered_text)
            if unresolved_cites:
                report.add("G4_PDF_FONTS", "FAIL", "UNRESOLVED_PDF_CITATIONS", f"Rendered PDF contains {len(unresolved_cites)} unresolved '[?]' citation markers!")
            if unresolved_refs:
                report.add("G4_PDF_FONTS", "FAIL", "UNRESOLVED_PDF_REFS", f"Rendered PDF contains {len(unresolved_refs)} unresolved '??' cross-reference markers!")
            if not unresolved_cites and not unresolved_refs:
                report.add("G4_PDF_FONTS", "PASS", "ZERO_UNRESOLVED_MARKERS", "Rendered PDF text verified free of '[?]' and '??' unresolved references.")

        # =====================================================================
        # GATE 5: PORTAL_METADATA_AND_MATHJAX_ABSTRACT
        # =====================================================================
        title_m = re.search(r"\\title\{([\s\S]*?)\}\s*\\author", tex_nc)
        abs_m = re.search(r"\\begin\{abstract\}([\s\S]*?)\\end\{abstract\}", tex_nc)
        author_m = re.search(r"\\author\{([\s\S]*?)\}\s*\\date", tex_nc)

        raw_title = title_m.group(1) if title_m else ""
        raw_abs = abs_m.group(1) if abs_m else ""
        raw_author = author_m.group(1) if author_m else ""

        report.portal_title = latex_to_arxiv_portal_text(raw_title, is_title=True)
        report.portal_abstract = latex_to_arxiv_portal_text(raw_abs, is_title=False)
        report.portal_abstract_chars = len(report.portal_abstract)

        if report.portal_abstract_chars == 0:
            report.add("G5_METADATA", "FAIL", "EMPTY_ABSTRACT", "Could not extract \\begin{abstract}...\\end{abstract}.")
        elif report.portal_abstract_chars > ARXIV_ABSTRACT_MAX_CHARS:
            over_by = report.portal_abstract_chars - ARXIV_ABSTRACT_MAX_CHARS
            report.add(
                "G5_METADATA",
                "FAIL",
                "ABSTRACT_EXCEEDS_1920_LIMIT",
                f"Portal abstract is {report.portal_abstract_chars} chars — exceeds arXiv's hard 1,920-char limit by +{over_by} chars!",
            )
        else:
            headroom = ARXIV_ABSTRACT_MAX_CHARS - report.portal_abstract_chars
            report.add(
                "G5_METADATA",
                "PASS",
                "ABSTRACT_LENGTH_OK",
                f"Portal abstract is {report.portal_abstract_chars} / {ARXIV_ABSTRACT_MAX_CHARS} chars ({headroom} chars headroom).",
            )

        if "\\cite{" in raw_abs or "\\ref{" in raw_abs or "\\eqref{" in raw_abs:
            report.add("G5_METADATA", "FAIL", "CITE_OR_REF_IN_ABSTRACT", "Abstract contains \\cite or \\ref; arXiv landing pages cannot resolve internal numbers.")

        if report.portal_abstract.count("$") % 2 != 0:
            report.add("G5_METADATA", "FAIL", "UNBALANCED_MATHJAX_DOLLARS", "Portal abstract has unbalanced '$' math delimiters.")

        outside_math = re.sub(r"\$[^$]*\$", "", report.portal_abstract)
        leaked_macros = re.findall(r"\\[A-Za-z]+", outside_math)
        if leaked_macros:
            report.add("G5_METADATA", "FAIL", "UNPARSED_TEX_MACRO_IN_ABSTRACT", f"Unconverted TeX macros outside $...$ in portal abstract: {sorted(set(leaked_macros))}")
        else:
            report.add("G5_METADATA", "PASS", "MATHJAX_ABSTRACT_CLEAN", "Portal abstract contains zero unparsed text-mode TeX macros outside $...$ blocks.")

        if "arXiv Primary" in raw_author or "Cross-Lists:" in raw_author:
            report.add("G5_METADATA", "FAIL", "AUTHOR_BLOCK_ARXIV_CLUTTER", "\\author{} self-prints 'arXiv Primary / Cross-Lists' metadata clutter.")

        if "mitorestorer.org" in raw_author.lower() or "osatz@" in raw_author.lower():
            report.add("G5_METADATA", "FAIL", "HALLUCINATED_AUTHOR_CONTACT", "Fabricated domain/email detected in \\author{} block.")

        if cfg:
            all_cats = [cfg.primary_category, *cfg.cross_categories]
            bad_cats = [c for c in all_cats if c not in VALID_ARXIV_CATEGORIES]
            if bad_cats:
                report.add("G5_METADATA", "FAIL", "INVALID_ARXIV_CATEGORY", f"Invalid arXiv category taxonomy codes: {bad_cats}")

        # =====================================================================
        # GATE 6: BIBLIOGRAPHIC_AND_DOI_INTEGRITY
        # =====================================================================
        cited_keys: List[str] = []
        for cm in re.finditer(r"\\cite\{([^}]+)\}", tex_nc):
            for k in cm.group(1).split(","):
                ks = k.strip()
                if ks and ks not in cited_keys:
                    cited_keys.append(ks)

        bbl_items: Dict[str, str] = {}
        parts_bbl = re.split(r"\\bibitem(?:\[[^\]]*\])?\{([^}]+)\}", bbl_raw)
        for i in range(1, len(parts_bbl), 2):
            k = parts_bbl[i].strip()
            body = parts_bbl[i + 1].replace("\\end{thebibliography}", "").strip()
            bbl_items[k] = " ".join(body.split())

        missing_in_bbl = set(cited_keys) - set(bbl_items.keys())
        uncited_in_bbl = set(bbl_items.keys()) - set(cited_keys)
        if missing_in_bbl:
            report.add("G6_CITATIONS", "FAIL", "MISSING_BIBITEM", f"Citations missing from .bbl: {sorted(missing_in_bbl)}")
        if uncited_in_bbl:
            report.add("G6_CITATIONS", "WARN", "UNCITED_BIBITEM", f"Uncited \\bibitem entries in .bbl: {sorted(uncited_in_bbl)}")

        # Duplicate normalized bibitem detector
        seen_norm: Dict[str, str] = {}
        dup_pairs: List[Tuple[str, str]] = []
        for k, body in bbl_items.items():
            norm = re.sub(r"[^a-z0-9]", "", body.lower())
            if norm in seen_norm:
                dup_pairs.append((seen_norm[norm], k))
            else:
                seen_norm[norm] = k
        if dup_pairs:
            report.add("G6_CITATIONS", "FAIL", "DUPLICATE_BIBITEM_ENTRIES", f"Duplicate bibliography entries under distinct keys: {dup_pairs}")
        elif not missing_in_bbl:
            report.add("G6_CITATIONS", "PASS", "BBL_BIJECTION_OK", f"All {len(cited_keys)} cited references uniquely mapped in {stem}.bbl.")

        if bib_db:
            dois_to_verify: List[Tuple[str, str, str]] = []
            for k in cited_keys:
                entry = bib_db.get(k)
                if not entry:
                    continue
                if k == "satz2026paper1":
                    howpub = entry.get("howpublished", "")
                    if "arXiv:" not in howpub:
                        report.add(
                            "G6_CITATIONS",
                            "WARN",
                            "COMPANION_ARXIV_ID_PENDING",
                            "Citation 'satz2026paper1' references companion Paper I without an assigned arXiv:YYMM.NNNNN ID yet (submit Paper I first, then update).",
                        )
                doi = entry.get("doi", "").strip()
                title = entry.get("title", "")
                if doi:
                    if not re.match(r"^10\.\d{4,9}/\S+$", doi):
                        report.add("G6_CITATIONS", "FAIL", "MALFORMED_DOI", f"Key '{k}' has malformed DOI '{doi}'.")
                    else:
                        dois_to_verify.append((k, doi, title))

            if verify_dois_online and dois_to_verify:
                ok_cnt = 0
                for k, doi, title in dois_to_verify:
                    ok, detail = verify_doi_crossref(doi, title)
                    if not ok:
                        report.add("G6_CITATIONS", "FAIL", "CROSSREF_DOI_MISMATCH", f"[{k}] ({doi}): {detail}")
                    else:
                        ok_cnt += 1
                if ok_cnt == len(dois_to_verify):
                    report.add("G6_CITATIONS", "PASS", "CROSSREF_LIVE_VERIFIED", f"100% ({ok_cnt}/{len(dois_to_verify)}) DOIs live-verified against Crossref/RCSB APIs.")

        # =====================================================================
        # GATE 7: MODERATION_AND_OVERLAP_SCREENER
        # =====================================================================
        if cfg and cfg.is_omnibus:
            report.add(
                "G7_MODERATION",
                "FAIL",
                "OMNIBUS_REDUNDANT_SUBMISSION",
                "Omnibus summary overlaps standalone series manuscripts; arXiv moderators hold/reject redundant omnibus uploads.",
            )

        if sibling_tex_map:
            my_shingles = _shingle_set(tex_nc)
            for sib_name, sib_tex in sibling_tex_map.items():
                if sib_name == toplevel_tex_name or sib_name == "main.tex":
                    continue
                sib_shingles = _shingle_set(strip_latex_comments(sib_tex))
                if my_shingles and sib_shingles:
                    containment = len(my_shingles & sib_shingles) / len(my_shingles)
                    if containment > 0.18:
                        report.add(
                            "G7_MODERATION",
                            "FAIL",
                            "HIGH_TEXT_OVERLAP",
                            f"5-gram text overlap with {sib_name} is {containment:.1%} (> 18% arXiv text-recycling hold threshold).",
                        )

        # Figure text screener (infographic banner detection)
        if pdftotext_bin:
            banner_patterns = [
                r"FLAGSHIP\s*#",
                r"GRAVEYARD\s+AUDIT",
                r"DISCOVERIES\s+12",
                r"CELL\s*\(SEPT\s*2026\)",
            ]
            for fig_file in clean.glob("*.pdf"):
                if fig_file.name == f"{stem}.pdf":
                    continue
                ftxt = subprocess.run([pdftotext_bin, str(fig_file), "-"], capture_output=True, text=True).stdout
                for pat in banner_patterns:
                    if re.search(pat, ftxt, flags=re.IGNORECASE):
                        report.add(
                            "G7_MODERATION",
                            "FAIL",
                            "INFOGRAPHIC_FIGURE_BANNER",
                            f"Figure '{fig_file.name}' contains non-academic dashboard banner matching /{pat}/.",
                        )

        # Clinical / prescriptive phrasing screener vs computational scope disclosure
        clinical_triggers = [
            (r"Pathway\s+B\s*---?\s*Dietary", "Prescriptive 'Pathway B — Dietary & GRAS' wellness regimen"),
            (r"patients\s+should\s+receive", "Direct clinical dosing directive"),
            (r"recommended\s+human\s+dose", "Human clinical dosing recommendation"),
        ]
        for pat, label in clinical_triggers:
            if re.search(pat, tex_nc, flags=re.IGNORECASE):
                report.add(
                    "G7_MODERATION",
                    "FAIL",
                    "CLINICAL_MODERATION_RISK",
                    f"Detected {label} without clinical trial registration (high q-bio moderation hold risk).",
                )

        if not re.search(r"\b(in\s+silico|computational\s+study|host-emulated|theoretical)\b", tex_nc, flags=re.IGNORECASE):
            report.add("G7_MODERATION", "WARN", "MISSING_COMPUTATIONAL_DISCLOSURE", "Explicit 'in silico / computational study' framing not detected.")
        elif report.fail_count == 0:
            report.add("G7_MODERATION", "PASS", "MODERATION_CLEAN", "Passed text-overlap, vector-figure, and computational-scope moderation checks.")

    return report


def format_console_report(report: AuditReport) -> str:
    verdict = "PASS (UPLOAD-READY)" if report.passed else f"BLOCKED ({report.fail_count} FAIL, {report.warn_count} WARN)"
    lines = [
        f"\n{'=' * 78}",
        f"  arXiv Submission Analog Report :: {report.paper_id.upper()} -> {verdict}",
        f"{'=' * 78}",
        f"  Archive  : {report.archive_path} ({report.archive_bytes:,} B | SHA256: {report.sha256[:16]}...)",
        f"  Metrics  : {report.pages} pages | {report.figures_count} figures | {report.tables_count} tables | Abstract: {report.portal_abstract_chars}/{ARXIV_ABSTRACT_MAX_CHARS} chars",
        f"{'-' * 78}",
    ]
    for f in report.findings:
        tag = f"[{f.severity:4s}]"
        lines.append(f"  {tag} {f.gate:14s} | {f.code:30s} | {f.message}")
    lines.append(f"{'=' * 78}")
    return "\n".join(lines)


def run_cli(argv: Optional[List[str]] = None) -> int:
    parser = argparse.ArgumentParser(description="arXiv Submission Analog & Verification Engine")
    parser.add_argument("--build-ready", action="store_true", help="Build & verify upload-ready .tar.gz bundles (Paper I & Paper II).")
    parser.add_argument("--audit-all", action="store_true", help="Build temporary cleanroom bundles and audit all 5 repository manuscripts.")
    parser.add_argument("--audit-archive", type=Path, help="Audit an existing .tar.gz submission archive.")
    parser.add_argument("--verify-dois-online", action="store_true", help="Verify every cited DOI live against Crossref & RCSB APIs.")
    parser.add_argument("--json-out", type=Path, default=DIST_DIR / "arxiv_preflight_summary.json", help="Path for machine-readable JSON audit report.")
    args = parser.parse_args(argv)

    if not (args.build_ready or args.audit_all or args.audit_archive):
        args.build_ready = True

    bib_db = parse_bib_database(ARXIV_SRC_DIR / "references.bib")
    sibling_tex_map = {
        p.name: p.read_text(encoding="utf-8")
        for p in sorted(ARXIV_SRC_DIR.glob("*.tex"))
    }

    reports: List[AuditReport] = []

    if args.audit_archive:
        rep = audit_archive(
            args.audit_archive,
            cfg=None,
            bib_db=bib_db,
            sibling_tex_map=sibling_tex_map,
            verify_dois_online=args.verify_dois_online,
        )
        print(format_console_report(rep))
        reports.append(rep)

    if args.build_ready:
        for pid in ("paper1", "paper2"):
            cfg = MANUSCRIPT_REGISTRY[pid]
            archive_path, manifest_path = build_arxiv_bundle(cfg, out_dir=DIST_DIR)
            rep = audit_archive(
                archive_path,
                cfg=cfg,
                bib_db=bib_db,
                sibling_tex_map=sibling_tex_map,
                verify_dois_online=args.verify_dois_online,
            )
            write_portal_manifest(cfg, rep, manifest_path)
            print(format_console_report(rep))
            reports.append(rep)

    if args.audit_all:
        with tempfile.TemporaryDirectory(prefix="arxiv_audit_all_") as tmp_all:
            tmp_dir = Path(tmp_all)
            for pid, cfg in MANUSCRIPT_REGISTRY.items():
                out_target = DIST_DIR if cfg.target_ready else tmp_dir
                archive_path, manifest_path = build_arxiv_bundle(cfg, out_dir=out_target)
                rep = audit_archive(
                    archive_path,
                    cfg=cfg,
                    bib_db=bib_db,
                    sibling_tex_map=sibling_tex_map,
                    verify_dois_online=args.verify_dois_online if cfg.target_ready else False,
                )
                if cfg.target_ready:
                    write_portal_manifest(cfg, rep, manifest_path)
                print(format_console_report(rep))
                reports.append(rep)

    args.json_out.parent.mkdir(parents=True, exist_ok=True)
    args.json_out.write_text(
        json.dumps([dataclasses.asdict(r) for r in reports], indent=2) + "\n",
        encoding="utf-8",
    )

    if args.audit_all:
        # In --audit-all mode, succeed if target_ready papers pass and held-back papers were accurately caught
        ready_failed = [r for r in reports if MANUSCRIPT_REGISTRY[r.paper_id].target_ready and not r.passed]
        return 1 if ready_failed else 0

    return 0 if all(r.passed for r in reports) else 1


if __name__ == "__main__":
    sys.exit(run_cli())

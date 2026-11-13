# MitoRestorer Segregated arXiv Submission Guide & Cleanroom Verification Engine

To comply with arXiv Submission 1.5+ (`AutoTeX` / TeX Live `pdflatex`) packaging, font-embedding, abstract-length (`<= 1,920` chars), citation-authenticity, and category-moderation policies, the `MitoRestorer` program separates its manuscripts into standalone uploads verified by the automated **arXiv Submission Analog** ([`tools/arxiv_submission_analog.py`](file:///Users/o2satz/antigravity/mitorestorer/tools/arxiv_submission_analog.py)).

---

## 1. Upload-Ready Bundles (`dist/`) & Staged Submission Schedule

| Priority & Status | Upload Archive & Portal Manifest | Manuscript & Verified v2 Scope | Primary Category | Cross-Lists | Cleanroom Analog Verdict |
| :--- | :--- | :--- | :--- | :--- | :--- |
| **1st — UPLOAD TODAY** | **Archive**: [`dist/arxiv_paper1.tar.gz`](file:///Users/o2satz/antigravity/mitorestorer/dist/arxiv_paper1.tar.gz)<br>**Metadata**: [`dist/arxiv_paper1_portal_metadata.txt`](file:///Users/o2satz/antigravity/mitorestorer/dist/arxiv_paper1_portal_metadata.txt)<br>**PDF**: [`paper1_core_control_and_coral_npu.pdf`](file:///Users/o2satz/antigravity/mitorestorer/docs/arxiv_submission/paper1_core_control_and_coral_npu.pdf) | **Paper I**: *Nondimensional Staged Feedback Control of Mitochondrial Heteroplasmy Dynamics with Host-Emulated Verification on an Ultra-Low-Power `RV32IMF_Zve32x` Neural Processing Unit*<br>• Nagumo safe-set forward invariance ($\mathrm{EP} \le 1.0456 < 1.150$, $\lambda_{\min} = 0.01085\text{ d}^{-1}$)<br>• Exact selection collapse ($0.0000$ error) & non-homothetic scaling ($R^2 = 0.9828$, error $0.0493 \le 0.1216$)<br>• 300-var discrete-adjoint PMP ($J_{\mathrm{PMP}}^* = 0.7262$ vs. $J_{\mathrm{thresh}} = 0.8875$), Paulin Markov persistence bound, `3,696 B` L1 DTCM SRAM, `Q16.16/Q8.24` $\Sigma\Delta$ RK4 (`0.0534%–0.0563%` error) | **`eess.SY`** | `q-bio.QM`, `math.OC`, `cs.AR` | **`PASS`**<br>(`0 FAIL, 0 WARN`; `19/19` DOIs live-verified on Crossref; `1,873/1,920` abstract chars) |
| **2nd — Upload After Paper I Gets `arXiv:YYMM.NNNNN`** | **Archive**: [`dist/arxiv_paper2.tar.gz`](file:///Users/o2satz/antigravity/mitorestorer/dist/arxiv_paper2.tar.gz)<br>**Metadata**: [`dist/arxiv_paper2_portal_metadata.txt`](file:///Users/o2satz/antigravity/mitorestorer/dist/arxiv_paper2_portal_metadata.txt)<br>**PDF**: [`paper2_ex_vivo_cart_bioreactor.pdf`](file:///Users/o2satz/antigravity/mitorestorer/docs/arxiv_submission/paper2_ex_vivo_cart_bioreactor.pdf) | **Paper II**: *In Silico Closed-Loop Metabolic Gating of Glycolytic-to-Oxidative Transition in Ex Vivo CAR-T Perfusion Bioreactors Guided by Structural Bioinformatics and Optical Redox Observability*<br>• Atomic `PDK1` (`2Q8G`, `Ser75` `2.50 Å`) & `SIRT1` (`4ZZJ`, `Glu230` `3.47 Å`) + AlphaFold v6 `FOXO1` (`pLDDT 50.96`) & `PGC-1α` (`pLDDT 52.74`)<br>• 168 h Bioreactor-BCR (`62.80x` fold, `73.20%` $\mathrm{T}_{\mathrm{SCM}}$, `SRC = 4.87x`, `h = 6.73%`)<br>• $N = 1{,}000$ oncology donor cohort: cell-density-invariant `ORR` lifts 168 h pass rate `40.5% -> 58.7%` ($p = 8.04 \times 10^{-4}$ vs. 68 h timer), `59.4%` under biofouling ($p = 5.22 \times 10^{-129}$), `74.0%` adaptive harvest; `2,704 B` FPU-free C99 firmware | **`q-bio.CB`** | `q-bio.QM`, `eess.SY` | **`PASS`**<br>(`0 FAIL, 1 WARN` pending Paper I `arXiv` ID in `satz2026paper1`; `17/17` DOIs live-verified; `1,907/1,920` chars) |
| **3rd — Hold Pending v2 Benchmark & Figure Overhaul** | [`paper3_apis_pollinator_smart_hive.tex`](file:///Users/o2satz/antigravity/mitorestorer/docs/arxiv_submission/paper3_apis_pollinator_smart_hive.tex) | **Paper III**: *Apis mellifera* Flight-Muscle Holobiont Controller (`Apis-BCR`) | `q-bio.QM` | `q-bio.PE`, `eess.SY` | **Held** (needs `generate_paper3_verified_benchmarks.py` + Type-42 `matplotlib` figures) |
| **4th — Hold Pending v2 PDE/QSP & Clinical De-Risking** | [`paper4_neuro_qsp_and_axonal_pde.tex`](file:///Users/o2satz/antigravity/mitorestorer/docs/arxiv_submission/paper4_neuro_qsp_and_axonal_pde.tex) | **Paper IV**: Post-Mitotic Nigrostriatal QSP & Axonal Reaction-Diffusion PDE (`Neuro-BCR`) | `q-bio.NC` | `q-bio.QM`, `q-bio.TO` | **Held** (needs Proposition 2 fix, `generate_paper4_verified_benchmarks.py`, Type-42 figures, non-prescriptive QSP framing) |
| **Repository Series Overview (Do Not Submit to arXiv)** | [`main.pdf`](file:///Users/o2satz/antigravity/mitorestorer/docs/arxiv_submission/main.pdf) | **Omnibus Whitepaper**: Retained on GitHub / Web Reader only (arXiv moderators reject overlapping omnibus uploads). | — | — | **Blocked by Design** (`OMNIBUS_REDUNDANT_SUBMISSION`) |

---

## 2. How to Use the arXiv Submission Analog (`tools/arxiv_submission_analog.py`)

[`tools/arxiv_submission_analog.py`](file:///Users/o2satz/antigravity/mitorestorer/tools/arxiv_submission_analog.py) replicates all seven verification gates of arXiv's Submission 1.5+ / `AutoTeX` and moderation pipeline inside an isolated `/tmp` cleanroom with zero access to parent directories or `.bib` files:

1. **`G1_ARCHIVE` (Archive Quarantine & `00README.json` v1.5 Schema)**: Enforces `< 50 MB` archive / `< 10 MB` file limits, blocks `..` traversal, symlinks, macOS `._*` junk, `.aux`/`.log` leaks, pre-compiled root `.pdf`s, and validates `00README.json` (`compiler: pdflatex`).
2. **`G2_TEX_STATIC` (TeX Live `pdflatex` Static Lint)**: Checks `\pdfoutput=1`, `utf8` `inputenc` + `T1` `fontenc`, blocks `--shell-escape` / `\write18` / XeTeX-only packages, and enforces 100% 7-bit ASCII across `.tex` and `.bbl`.
3. **`G3_AUTOTEX` (Hermetic Cleanroom Compilation)**: Compiles the unpacked `.tar.gz` in an isolated temp folder without `references.bib` (via `\input{<stem>.bbl}`) and checks Overfull `\hbox` margin tolerances.
4. **`G4_PDF_FONTS` (Poppler `pdfinfo` / `pdffonts` / `pdftotext` Forensics)**: Verifies every font object is embedded as crisp vector Type-1C or TrueType (`Type-42`), blocks `Type 3` bitmap fonts, and scans rendered PDF text for `[?]` or `??` unresolved markers.
5. **`G5_METADATA` (Step-3 Portal & MathJax Abstract Validator)**: Compiles LaTeX into arXiv's MathJax landing-page format, enforces the **hard `1,920`-character abstract cutoff**, balanced `$...$` math, zero `\cite`/`\ref`/unparsed text macros in the abstract, and clean `\author{}` blocks (`satzwa@temple.edu`).
6. **`G6_CITATIONS` (`.bbl` Bijection & Live Crossref/RCSB DOI Verification)**: Verifies 1-to-1 `\cite` $\leftrightarrow$ `\bibitem` mapping, catches duplicate bibliography entries, and (with `--verify-dois-online`) queries `api.crossref.org` and `data.rcsb.org` to verify live title/DOI match.
7. **`G7_MODERATION` (Shingle Overlap, Figure Banner & Clinical Claim Screener)**: Computes 5-gram shingle overlap across companion papers, scans figure PDFs for non-academic banners, and screens for unverified clinical dosing language.

### Commands

```bash
# 1. Build upload-ready .tar.gz archives + Step-3 copy-paste manifests in dist/ and run offline 7-gate verification:
python3 tools/arxiv_submission_analog.py --build-ready

# 2. Run with live Crossref & RCSB API verification of every cited DOI:
python3 tools/arxiv_submission_analog.py --build-ready --verify-dois-online

# 3. Audit an arbitrary .tar.gz archive exactly as uploaded to arxiv.org/submit:
python3 tools/arxiv_submission_analog.py --audit-archive dist/arxiv_paper1.tar.gz

# 4. Run full readiness matrix across all 5 manuscripts in docs/arxiv_submission/:
python3 tools/arxiv_submission_analog.py --audit-all
```

---

## 3. Step-by-Step Upload Checklist for Paper I Today (`arxiv.org/submit`)

1. Sign in at `https://arxiv.org/user` using **`satzwa@temple.edu`**, confirm your ORCID iD is linked, and click **Start New Submission**.
2. **Step 1 (Category & License)**:
   * Primary category: **`eess.SY`** (Systems and Control)
   * Cross-list categories: **`q-bio.QM`**, **`math.OC`**, **`cs.AR`**
   * License: **`CC BY 4.0`**
3. **Step 2 (Upload Files)**:
   * Upload **`dist/arxiv_paper1.tar.gz`** (`168.9 KiB`, 5 files: `00README.json`, `.tex`, `.bbl`, and 3 vector `.pdf` figures).
   * Click **Process Files** and inspect arXiv's generated PDF preview (`15 pages`).
4. **Step 3 (Metadata)**:
   * Open [`dist/arxiv_paper1_portal_metadata.txt`](file:///Users/o2satz/antigravity/mitorestorer/dist/arxiv_paper1_portal_metadata.txt) and copy-paste **Title**, **Authors** (`Wayne A. Satz`), **Abstract** (`1,873 / 1,920` chars), **Comments**, **MSC-class**, and **ACM-class** directly into the portal fields.
5. **After Paper I Announces (`arXiv:2611.XXXXX`)**:
   * Update `satz2026paper1` in [`references.bib`](file:///Users/o2satz/antigravity/mitorestorer/docs/arxiv_submission/references.bib#L333-L338) with `arXiv:2611.XXXXX [eess.SY]`, re-run `python3 tools/arxiv_submission_analog.py --build-ready` (`0 FAIL, 0 WARN`), and upload **`dist/arxiv_paper2.tar.gz`** under **`q-bio.CB`**.

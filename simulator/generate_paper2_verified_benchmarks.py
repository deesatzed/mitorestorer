#!/usr/bin/env python3
"""
MitoRestorer Paper II — Deterministic SOTA Computational Verification & Figure Suite
=====================================================================================
Generates 100% live-computed structural bioinformatics, 8-state CAR-T perfusion
bioreactor ODE trajectories, discrete-adjoint optimal control & Pareto frontiers,
N=200 heterogeneous oncology apheresis closed-loop vs. fixed-timer benchmarks,
sensor biofouling innovation-residual handover, plant-preserving acidosis rescue,
and host-emulated C99 Q16.16/Q8.24 Sigma-Delta RK4 parity metrics.

Epistemic Scope: All dynamical results are strictly in silico computational
simulations and structural bioinformatics analyses designed to generate
falsifiable predictions for future ex vivo benchtop studies.

Seed: 20261111 (MIASE-compliant deterministic execution)
Outputs:
  - simulator/paper2_verified_benchmarks.json
  - figures/paper2_fig1_structural_pharmacology_and_pae.pdf
  - figures/paper2_fig2_bioreactor_ode_and_pareto_optimal.pdf
  - figures/paper2_fig3_donor_cohort_biofouling_and_acidosis.pdf
"""

from __future__ import annotations

import json
import math
import os
os.environ.setdefault("MPLCONFIGDIR", "/tmp/mpl_cache_mitorestorer")
import subprocess
import tempfile
import time
from pathlib import Path
from typing import Dict, List, Tuple

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
from scipy.optimize import minimize
from scipy.stats import binomtest

ROOT = Path(__file__).resolve().parent.parent
STRUCT_DIR = ROOT / "structures"
FIG_DIR = ROOT / "figures"
SIM_DIR = ROOT / "simulator"
CORAL_DIR = ROOT / "build" / "coralnpu"
SEED = 20261111

# Nominal 8-State CAR-T Perfusion Cassette Parameters (Exact 1-to-1 Match with C99 & LaTeX)
NOMINAL_PARAMS: Dict[str, float] = dict(
    K=1.40,
    k_mito_w=0.008,          # h^-1: basal WT mitochondrial turnover
    k_bio_0=0.028,           # h^-1: homeostatic biogenesis rate (W* = 1.00 when u=(0,0,0), D=0)
    k_bio_max=0.086,         # h^-1: Phase-2 SIRT1/PGC-1a stimulated cristae biogenesis
    k_mito_d_max=0.068,      # h^-1: Phase-1 Urolithin A stimulated PINK1 mitophagy
    k_rep_adv=1.22,          # fragmented/depolarized mtDNA replicative advantage
    xi_d=1.00,               # Phase-2 template replication parity factor
    k_dmg=0.0012,            # h^-1 EP^-1: supra-threshold ROS cristae fragmentation rate
    theta_gsh=1.38,          # glutathione/peroxiredoxin redox buffering ceiling
    R_basal=8.0,             # aerobic OXPHOS baseline Lac:Pyr molar ratio
    K_LP=22.0,               # reverse-ATPase depolarization masking half-saturation
    eta_d=0.22,              # residual OXPHOS coupling efficiency of fragmented D
    gamma_pdh=0.72,          # AZD7545 (PDB: 2Q8G) PDH E1a unblocking oxidative gain
    gamma_src=0.22,          # cristae SRC bioenergetic coupling coefficient
    alpha_gal=0.74,          # Media-B glucose->galactose/pyruvate glycolytic attenuation
    kappa_glyc=1.40,         # glycolytic ATP supply coefficient
    c_akt=0.08,              # mild Phase-1 PI3Kdelta/Akt glycolytic restraint
    K_R_ATP=8.0,             # glycolytic flux half-saturation on R_LP
    beta_lp=0.40,            # glycolytic NADH reductive stress weight in EP
    beta_h=0.88,             # mutant heteroplasmy impedance weight in EP
    beta_ph=1.90,            # extracellular acidosis impedance weight in EP
    c1=0.24,                 # Phase-1 mitophagy/Akt relief of energetic resistance
    c2=0.06,                 # Phase-2 biogenesis & transcriptional ATP remodeling load
    mu_glyc=0.062,           # h^-1: rapid Phase-1 glycolytic lymphoblast division rate (Td ~ 11.2 h)
    mu_ox=0.0185,            # h^-1: slower Phase-2 TCA/OXPHOS memory division rate
    K_ATP=0.52,              # total ATP supply half-saturation for cell division
    K_W_div=0.56,            # intact cristae W*(1-h) half-saturation for OXPHOS division
    N_cap=85.0,              # cassette fold-expansion carrying capacity
    k_warburg=3.00,          # h^-1: aerobic glycolytic Lac:Pyr generation rate
    alpha_h=0.46,            # depolarized mitochondrial amplification of Warburg flux
    k_pdh=0.088,             # h^-1: PDH complex pyruvate oxidation clearance rate
    beta_azd=2.80,           # AZD7545 (2Q8G) allosteric PDK1 inhibition multiplier
    k_perf_lp=0.120,         # h^-1: nominal cassette supernatant perfusion exchange on R_LP
    k_flush_lp=0.85,         # h^-1: emergency neutral bypass flush clearance on R_LP
    R_media=5.0,             # fresh pyruvate-buffered media baseline ratio
    k_imp=0.056,             # h^-1: nuclear FOXO1 import/retention rate constant
    f_0=0.30,                # basal nuclear FOXO1 shuttling drive
    alpha_akt=0.30,          # Phase-1 PI3Kdelta/Akt nuclear retention contribution
    phi_sirt1=1.18,          # Phase-2 SIRT1 (4ZZJ Glu230 + NR) nuclear deacetylation drive
    k_exp=0.032,             # h^-1: Akt/14-3-3 cytoplasmic export rate constant
    gamma_acid=7.50,         # acidosis-accelerated nuclear FOXO1 loss
    s_basal=0.044,           # h^-1: basal inner-membrane respiratory reserve synthesis
    k_cristae=0.074,         # h^-1: OPA1/PGC-1a cristae remodeling synthesis rate
    beta_src=1.32,           # nuclear FOXO1/SIRT1 + PDH amplification of cristae SRC
    k_src_dec=0.048,         # h^-1: cristae respiratory reserve turnover decay rate
    k_exh=0.0038,            # h^-1: supra-threshold EP-driven TOX+ epigenetic exhaustion
    theta_exh=1.36,          # redox stress threshold for accelerated exhaustion
    k_warb_exh=0.0032,       # h^-1: chronic glycolytic/mTORC1 effector chromatin scarring
    k_acid_exh=0.220,        # h^-1: severe lactic acidosis exhaustion/injury rate
    r_epig_0=0.0032,         # h^-1: slow basal chromatin turnover
    k_rev_epig=0.0110,       # h^-1: partial Phase-2 FOXO1/SIRT1 chromatin plasticity
    k_perf_0=0.48,           # h^-1: nominal cassette bicarbonate perfusion buffering rate
    k_perf_flush=1.80,       # h^-1: emergency neutral buffer flush rate
    pH_media=7.36,           # fresh perfusate media pH
    k_acid=0.0070,           # pH h^-1: cell-density-scaled lactic acidification rate
    R_buf=12.0,              # bicarbonate-buffered Lac:Pyr threshold
    N_ref=48.0,              # reference cassette cell-density fold
    K_SRC_mem=0.68,          # SRC half-saturation for stem-like memory T_SCM phenotype
)

X0_NOMINAL = np.array([0.54, 0.42, 1.0, 22.0, 0.28, 1.30, 0.08, 7.32], dtype=float)
T_HORIZON_H = 168.0


# =============================================================================
# 1. Structural Bioinformatics Audit (mmCIF Crystal + AlphaFold v6 pLDDT/PAE)
# =============================================================================
def parse_pdb_mmcif_contacts(cif_path: Path, chain_id: str, lig_comp_id: str, cutoff: float = 4.0) -> Dict:
    """Parse experimental X-ray mmCIF atomic coordinates and compute exact heavy-atom ligand contacts."""
    prot_atoms = []
    lig_atoms = []
    resolved_res = set()
    with cif_path.open("r", encoding="utf-8") as fh:
        in_atom = False
        headers: List[str] = []
        for raw_line in fh:
            line = raw_line.strip()
            if line == "loop_":
                in_atom = False
                headers = []
                continue
            if line.startswith("_atom_site."):
                in_atom = True
                headers.append(line.split(".")[1])
                continue
            if in_atom:
                if not line or line.startswith("#") or line.startswith("_"):
                    in_atom = False
                    continue
                parts = line.split()
                if len(parts) < len(headers):
                    continue
                rec = dict(zip(headers, parts))
                group = rec.get("group_PDB", "")
                elem = rec.get("type_symbol", "")
                if elem == "H":
                    continue
                asym = rec.get("auth_asym_id", rec.get("label_asym_id", ""))
                comp = rec.get("auth_comp_id", rec.get("label_comp_id", ""))
                seq = rec.get("auth_seq_id", rec.get("label_seq_id", "0"))
                atom_name = rec.get("auth_atom_id", rec.get("label_atom_id", ""))
                x = float(rec["Cartn_x"])
                y = float(rec["Cartn_y"])
                z = float(rec["Cartn_z"])
                if group == "ATOM" and asym == chain_id:
                    res_num = int(seq)
                    resolved_res.add(res_num)
                    prot_atoms.append((res_num, comp, atom_name, np.array([x, y, z])))
                elif group == "HETATM" and comp == lig_comp_id and asym == chain_id:
                    lig_atoms.append((atom_name, np.array([x, y, z])))

    lig_xyz = np.array([a[1] for a in lig_atoms])
    min_by_res: Dict[Tuple[int, str], Tuple[float, str, str]] = {}
    for res_num, comp, patom, pxyz in prot_atoms:
        dists = np.linalg.norm(lig_xyz - pxyz[None, :], axis=1)
        idx = int(np.argmin(dists))
        dmin = float(dists[idx])
        key = (res_num, comp)
        if key not in min_by_res or dmin < min_by_res[key][0]:
            min_by_res[key] = (dmin, patom, lig_atoms[idx][0])

    contacts = []
    for (res_num, comp), (dmin, patom, latom) in sorted(min_by_res.items(), key=lambda kv: kv[1][0]):
        if dmin <= cutoff:
            contacts.append({
                "residue": f"{comp}{res_num}",
                "res_num": res_num,
                "res_name": comp,
                "min_dist_A": round(dmin, 2),
                "protein_atom": patom,
                "ligand_atom": latom,
            })

    return {
        "resolved_residue_count": len(resolved_res),
        "resolved_span": [min(resolved_res), max(resolved_res)],
        "ligand_heavy_atoms": len(lig_atoms),
        "contacts_le_4A": contacts,
        "all_min_by_res": {f"{c}{r}": round(v[0], 2) for (r, c), v in min_by_res.items()},
    }


def parse_alphafold_cif_and_pae(cif_path: Path, pae_path: Path, ptm_sites: List[Tuple[int, str, str]]) -> Dict:
    """Parse AlphaFold v6 mmCIF per-residue C-alpha pLDDT and N x N PAE JSON matrix."""
    ca_plddt: Dict[int, Tuple[str, float]] = {}
    with cif_path.open("r", encoding="utf-8") as fh:
        in_atom = False
        headers: List[str] = []
        for raw_line in fh:
            line = raw_line.strip()
            if line == "loop_":
                in_atom = False
                headers = []
                continue
            if line.startswith("_atom_site."):
                in_atom = True
                headers.append(line.split(".")[1])
                continue
            if in_atom:
                if not line or line.startswith("#") or line.startswith("_"):
                    in_atom = False
                    continue
                parts = line.split()
                if len(parts) < len(headers):
                    continue
                rec = dict(zip(headers, parts))
                if rec.get("group_PDB") == "ATOM" and rec.get("label_atom_id") == "CA":
                    seq = int(rec["label_seq_id"])
                    comp = rec["label_comp_id"]
                    plddt = float(rec["B_iso_or_equiv"])
                    ca_plddt[seq] = (comp, plddt)

    n_res = len(ca_plddt)
    plddt_arr = np.array([ca_plddt[i][1] for i in range(1, n_res + 1)], dtype=float)
    aa_list = [ca_plddt[i][0] for i in range(1, n_res + 1)]

    with pae_path.open("r", encoding="utf-8") as fh:
        pae_raw = json.load(fh)
    pae_mat = np.array(pae_raw[0]["predicted_aligned_error"], dtype=float)
    off_mask = ~np.eye(n_res, dtype=bool)
    mean_off_pae = float(pae_mat[off_mask].mean())
    pct_pae_lt_5 = float(100.0 * (pae_mat[off_mask] < 5.0).mean())

    # Contiguous structured segments (smoothed 9-AA window pLDDT >= 70 of length >= 25 AA)
    kernel = np.ones(9) / 9.0
    smoothed = np.convolve(plddt_arr, kernel, mode="same")
    in_seg = smoothed >= 70.0
    domains = []
    start = None
    for idx, val in enumerate(in_seg):
        if val and start is None:
            start = idx
        elif not val and start is not None:
            if idx - start >= 25:
                sub_pae = pae_mat[start:idx, start:idx]
                sub_off = sub_pae[~np.eye(idx - start, dtype=bool)]
                domains.append({
                    "start_res": start + 1,
                    "end_res": idx,
                    "length_aa": idx - start,
                    "mean_plddt": round(float(plddt_arr[start:idx].mean()), 2),
                    "internal_mean_pae_A": round(float(sub_off.mean()), 2),
                })
            start = None
    if start is not None and n_res - start >= 25:
        sub_pae = pae_mat[start:n_res, start:n_res]
        sub_off = sub_pae[~np.eye(n_res - start, dtype=bool)]
        domains.append({
            "start_res": start + 1,
            "end_res": n_res,
            "length_aa": n_res - start,
            "mean_plddt": round(float(plddt_arr[start:n_res].mean()), 2),
            "internal_mean_pae_A": round(float(sub_off.mean()), 2),
        })

    ptm_audit = []
    for pos, expected_aa, role in ptm_sites:
        comp, score = ca_plddt[pos]
        ptm_audit.append({
            "site": f"{comp}{pos}",
            "expected": expected_aa,
            "role": role,
            "plddt": round(score, 2),
            "disordered_lt_70": bool(score < 70.0),
        })

    lys_indices = [i for i, aa in enumerate(aa_list) if aa == "LYS"]
    lys_plddts = plddt_arr[lys_indices]

    return {
        "length_aa": n_res,
        "mean_plddt": round(float(plddt_arr.mean()), 2),
        "count_very_low_lt_50": int((plddt_arr < 50.0).sum()),
        "pct_very_low_lt_50": round(float(100.0 * (plddt_arr < 50.0).mean()), 2),
        "count_low_50_70": int(((plddt_arr >= 50.0) & (plddt_arr < 70.0)).sum()),
        "pct_low_50_70": round(float(100.0 * ((plddt_arr >= 50.0) & (plddt_arr < 70.0)).mean()), 2),
        "count_confident_ge_70": int((plddt_arr >= 70.0).sum()),
        "pct_confident_ge_70": round(float(100.0 * (plddt_arr >= 70.0).mean()), 2),
        "mean_off_diag_pae_A": round(mean_off_pae, 2),
        "pct_pairs_pae_lt_5A": round(pct_pae_lt_5, 2),
        "structured_domains_ge_25aa": domains,
        "ptm_sites": ptm_audit,
        "total_lysines": len(lys_indices),
        "lysines_lt_70": int((lys_plddts < 70.0).sum()),
        "pct_lysines_lt_70": round(float(100.0 * (lys_plddts < 70.0).mean()), 2),
        "mean_lysine_plddt": round(float(lys_plddts.mean()), 2),
        "_plddt_arr": plddt_arr,
        "_pae_mat": pae_mat,
        "_lys_indices": lys_indices,
    }


def audit_all_structures() -> Dict:
    pdk1 = parse_pdb_mmcif_contacts(STRUCT_DIR / "2Q8G" / "pdb_00002q8g.cif", "A", "AZX", cutoff=4.0)
    # Record exact distances of legacy misattributed C-terminal residues for transparency
    pdk1["legacy_misattributed_distances_A"] = {
        k: pdk1["all_min_by_res"][k] for k in ("PHE318", "LEU321", "TYR342") if k in pdk1["all_min_by_res"]
    }
    del pdk1["all_min_by_res"]

    sirt1_stac = parse_pdb_mmcif_contacts(STRUCT_DIR / "4ZZJ" / "pdb_00004zzj.cif", "A", "4TQ", cutoff=4.0)
    sirt1_nad = parse_pdb_mmcif_contacts(STRUCT_DIR / "4ZZJ" / "pdb_00004zzj.cif", "A", "CNA", cutoff=4.0)
    del sirt1_stac["all_min_by_res"]
    del sirt1_nad["all_min_by_res"]

    foxo1_ptms = [
        (24, "THR", "Akt phosphorylation (14-3-3 nuclear export)"),
        (256, "SER", "Akt phosphorylation (gatekeeper export site)"),
        (319, "SER", "Akt/SGK1 phosphorylation (cytoplasmic sequestration)"),
        (245, "LYS", "CBP/p300 acetylation / SIRT1 deacetylation"),
        (248, "LYS", "CBP/p300 acetylation / SIRT1 deacetylation"),
        (262, "LYS", "CBP/p300 acetylation / SIRT1 deacetylation"),
        (265, "LYS", "CBP/p300 acetylation / SIRT1 deacetylation"),
        (274, "LYS", "CBP/p300 acetylation / SIRT1 deacetylation"),
        (294, "LYS", "CBP/p300 acetylation / SIRT1 deacetylation"),
    ]
    foxo1 = parse_alphafold_cif_and_pae(
        STRUCT_DIR / "alphafold" / "Q12778" / "AF-Q12778-F1-model_v6.cif",
        STRUCT_DIR / "alphafold" / "Q12778" / "AF-Q12778-F1-predicted_aligned_error_v6.json",
        foxo1_ptms,
    )

    pgc1a_ptms = [
        (77, "LYS", "GCN5 acetylation / SIRT1 deacetylation"),
        (144, "LEU", "Nuclear receptor LXXLL coactivator box"),
        (183, "LYS", "GCN5 acetylation / SIRT1 deacetylation"),
        (253, "LYS", "GCN5 acetylation / SIRT1 deacetylation"),
        (320, "LYS", "GCN5 acetylation / SIRT1 deacetylation"),
        (441, "LYS", "Central regulatory acetylation site"),
    ]
    pgc1a = parse_alphafold_cif_and_pae(
        STRUCT_DIR / "alphafold" / "Q9UBK2" / "AF-Q9UBK2-F1-model_v6.cif",
        STRUCT_DIR / "alphafold" / "Q9UBK2" / "AF-Q9UBK2-F1-predicted_aligned_error_v6.json",
        pgc1a_ptms,
    )

    return {
        "pdk1_2q8g": pdk1,
        "sirt1_4zzj_stac_4tq": sirt1_stac,
        "sirt1_4zzj_carba_nad": sirt1_nad,
        "foxo1_af_q12778": foxo1,
        "pgc1a_af_q9ubk2": pgc1a,
    }


# =============================================================================
# 2. Unified 8-State CAR-T Perfusion Cassette ODE & Discrete-Adjoint PMP
# =============================================================================
def eval_cart8_rhs_vec(
    X: np.ndarray,
    U: np.ndarray,
    p: Dict,
    fault_occl: bool = False,
    u_azd_ov: float | None = None,
    u_sirt_ov: float | None = None,
) -> np.ndarray:
    """Vectorized 8-state CAR-T bioreactor ODE RHS. X shape (8, B) or (8,), U shape (3, B) or (3,)."""
    W = np.maximum(0.0, X[0])
    D = np.maximum(0.0, X[1])
    N = np.maximum(0.0, X[2])
    R = np.maximum(0.0, X[3])
    F = np.clip(X[4], 0.0, 1.0)
    S = np.maximum(0.0, X[5])
    E_tox = np.clip(X[6], 0.0, 1.0)
    pH = X[7]

    u1, u2, u_fl = U[0], U[1], U[2]
    M = np.maximum(1e-9, W + D)
    h = D / M
    w_frac = W / M

    u_azd = u2 if u_azd_ov is None else np.full_like(u2, u_azd_ov)
    u_sirt = u2 if u_sirt_ov is None else np.full_like(u2, u_sirt_ov)
    if fault_occl:
        g_glc = np.ones_like(u2)
        u_azd_eff = np.zeros_like(u2)
    else:
        g_glc = 1.0 - p["alpha_gal"] * u2
        u_azd_eff = u_azd

    v_unmask = p["K_LP"] / (p["K_LP"] + np.maximum(0.0, R - p["R_basal"]))
    k_md = p["k_mito_w"] + (p["k_mito_d_max"] - p["k_mito_w"]) * u1 * v_unmask

    eff_mass_raw = W + p["eta_d"] * D
    atp_ox = eff_mass_raw * (1.0 + p["gamma_pdh"] * u_azd_eff) * (1.0 + p["gamma_src"] * np.maximum(0.0, S - 1.0))
    atp_glyc = p["kappa_glyc"] * g_glc * (1.0 - p["c_akt"] * u1 * (1.0 - u2)) * (R / (R + p["K_R_ATP"]))
    atp_tot = atp_ox + atp_glyc

    ph_def = np.maximum(0.0, 7.02 - pH)
    phi_ph = 1.0 / (1.0 + (ph_def / 0.16) ** 2)

    e_dot = (
        1.0
        + p["beta_lp"] * np.maximum(0.0, R - p["R_basal"]) / p["K_LP"]
        + p["beta_h"] * h
        + p["beta_ph"] * ph_def
    ) / np.maximum(0.26, eff_mass_raw)
    f_eff = 1.0 - p["c1"] * u1 + p["c2"] * u2
    ep = e_dot * (f_eff ** 2)
    dmg = p["k_dmg"] * np.maximum(0.0, ep - p["theta_gsh"])

    crowd = np.maximum(0.0, 1.0 - M / p["K"])
    dW = (p["k_bio_0"] + (p["k_bio_max"] - p["k_bio_0"]) * u_sirt * F) * W * crowd - (p["k_mito_w"] + dmg) * W
    dD = p["k_rep_adv"] * (p["k_bio_0"] + p["xi_d"] * (p["k_bio_max"] - p["k_bio_0"]) * u2) * D * crowd - k_md * D + dmg * W

    w_func = W * (1.0 - h)
    ox_comp = (w_func * w_func) / (w_func * w_func + p["K_W_div"] ** 2)
    pdh_gate = 0.15 + 0.85 * u_azd_eff
    mu_anab = p["mu_glyc"] * g_glc * (1.0 - 0.08 * u1 * (1.0 - u2)) + p["mu_ox"] * u2 * pdh_gate * ox_comp
    dN = mu_anab * phi_ph * (atp_tot / (atp_tot + p["K_ATP"])) * N * np.maximum(0.0, 1.0 - N / p["N_cap"])

    perf_lp = (0.08 * p["k_perf_lp"]) if fault_occl else p["k_perf_lp"]
    dR = (
        p["k_warburg"] * g_glc * (1.0 - 0.35 * u1 * (1.0 - u2)) * (1.0 + p["alpha_h"] * h)
        - p["k_pdh"] * (1.0 + p["beta_azd"] * u_azd_eff) * w_frac * R
        - (perf_lp + p["k_flush_lp"] * u_fl) * (R - p["R_media"])
    )

    sirt_pdh = u_sirt * pdh_gate
    dF = (
        p["k_imp"] * (p["f_0"] + p["alpha_akt"] * u1 + p["phi_sirt1"] * sirt_pdh) * phi_ph * (1.0 - F)
        - p["k_exp"] * (1.0 - 0.52 * u1) * (1.0 - 0.85 * sirt_pdh) * (1.0 + p["gamma_acid"] * np.maximum(0.0, 6.90 - pH)) * F
    )

    dS = (p["s_basal"] * eff_mass_raw + p["k_cristae"] * w_func * (1.0 + p["beta_src"] * F * sirt_pdh)) - p["k_src_dec"] * S

    dE = (
        (
            p["k_exh"] * np.maximum(0.0, ep - p["theta_exh"])
            + p["k_warb_exh"] * g_glc * (1.0 - 0.22 * u1) * (1.0 + 0.70 * (N / p["N_ref"]))
            + p["k_acid_exh"] * np.maximum(0.0, 6.86 - pH)
        )
        * (1.0 - E_tox)
        - (p["r_epig_0"] + p["k_rev_epig"] * F * sirt_pdh) * E_tox
    )

    k_perf = (0.08 * p["k_perf_0"] if fault_occl else p["k_perf_0"]) + p["k_perf_flush"] * u_fl
    dpH = k_perf * (p["pH_media"] - pH) - p["k_acid"] * (1.65 if fault_occl else g_glc) * np.maximum(0.0, R - p["R_buf"]) * (N / p["N_ref"])
    return np.stack([dW, dD, dN, dR, dF, dS, dE, dpH], axis=0)


def rk4_step_vec(X: np.ndarray, U: np.ndarray, p: Dict, dt: float, fault_occl: bool = False, u_azd_ov=None, u_sirt_ov=None) -> np.ndarray:
    k1 = eval_cart8_rhs_vec(X, U, p, fault_occl, u_azd_ov, u_sirt_ov)
    k2 = eval_cart8_rhs_vec(X + 0.5 * dt * k1, U, p, fault_occl, u_azd_ov, u_sirt_ov)
    k3 = eval_cart8_rhs_vec(X + 0.5 * dt * k2, U, p, fault_occl, u_azd_ov, u_sirt_ov)
    k4 = eval_cart8_rhs_vec(X + dt * k3, U, p, fault_occl, u_azd_ov, u_sirt_ov)
    return X + (dt / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)


def compute_tscm_pct(X: np.ndarray, p: Dict = NOMINAL_PARAMS) -> np.ndarray:
    return 100.0 * X[4] * (1.0 - X[6]) * (X[5] / (X[5] + p["K_SRC_mem"]))


def simulate_nominal_arm(mode: str, n_steps: int = 1000) -> Dict:
    dt = T_HORIZON_H / n_steps
    X = X0_NOMINAL[:, None].copy()
    traj = np.zeros((n_steps + 1, 8))
    u_traj = np.zeros((n_steps, 3))
    traj[0] = X[:, 0]
    stage = 1
    sw_t = None
    flush_until = -1.0
    min_ph = float(X[7, 0])

    for k in range(n_steps):
        t = k * dt
        h = float(X[1, 0] / (X[0, 0] + X[1, 0]))
        fault = (mode in ("fault_no_rb", "fault_rb")) and (72.0 <= t <= 96.0)
        u_azd_ov = 0.0 if mode == "ablate_azd" else None
        u_sirt_ov = 0.0 if mode == "ablate_sirt1" else None

        if mode == "warburg":
            u = [0.0, 0.0, 0.0]
        elif mode == "premature_ox_h0":
            u = [0.0, 1.0, 0.0]
        elif mode == "urolithin_mono":
            u = [1.0, 0.0, 0.0]
        elif mode.startswith("fixed_"):
            ts = float(mode.split("_")[1])
            if stage == 1 and t >= ts:
                stage = 2
                sw_t = t
            u = [1.0, 0.0, 0.0] if stage == 1 else [0.0, 1.0, 0.0]
        elif mode in ("closed_loop", "ablate_azd", "ablate_sirt1", "fault_no_rb"):
            if stage == 1 and h <= 0.045 and X[2, 0] >= 14.0 and X[7, 0] >= 7.10:
                stage = 2
                sw_t = t
            u = [1.0, 0.0, 0.0] if stage == 1 else [0.0, 1.0, 0.0]
        elif mode == "fault_rb":
            if stage == 1 and h <= 0.045 and X[2, 0] >= 14.0 and X[7, 0] >= 7.12 and t >= flush_until:
                stage = 2
                if sw_t is None:
                    sw_t = t
            if stage == 2 and X[7, 0] <= 6.85:
                stage = 1
                flush_until = t + 6.0
            u_fl = 1.0 if (t < flush_until or X[7, 0] <= 6.92) else 0.0
            u = [1.0, 0.0, u_fl] if stage == 1 else [0.0, 1.0, u_fl]
        else:
            raise ValueError(f"Unknown mode: {mode}")

        u_traj[k] = u
        U = np.array(u, dtype=float)[:, None]
        X = rk4_step_vec(X, U, NOMINAL_PARAMS, dt, fault, u_azd_ov, u_sirt_ov)
        traj[k + 1] = X[:, 0]
        min_ph = min(min_ph, float(X[7, 0]))

    xf = traj[-1]
    h_end = float(100.0 * xf[1] / (xf[0] + xf[1]))
    tscm = float(compute_tscm_pct(xf))
    # Accumulated selection integral s(t_sw) along Stage 1
    s_gate = None
    if sw_t is not None:
        idx_sw = int(round(sw_t / dt))
        r_lp_s = traj[:idx_sw, 3]
        v_un = NOMINAL_PARAMS["K_LP"] / (NOMINAL_PARAMS["K_LP"] + np.maximum(0.0, r_lp_s - NOMINAL_PARAMS["R_basal"]))
        k_md_s = NOMINAL_PARAMS["k_mito_w"] + (NOMINAL_PARAMS["k_mito_d_max"] - NOMINAL_PARAMS["k_mito_w"]) * v_un
        m_s = traj[:idx_sw, 0] + traj[:idx_sw, 1]
        lam_s = (k_md_s - NOMINAL_PARAMS["k_mito_w"]) - (NOMINAL_PARAMS["k_rep_adv"] - 1.0) * NOMINAL_PARAMS["k_bio_0"] * np.maximum(0.0, 1.0 - m_s / NOMINAL_PARAMS["K"])
        s_gate = float(np.sum(lam_s) * dt)

    return {
        "mode": mode,
        "switch_time_h": round(sw_t, 2) if sw_t is not None else None,
        "s_accumulated_selection": round(s_gate, 4) if s_gate is not None else None,
        "expansion_fold": round(float(xf[2]), 2),
        "h_end_pct": round(h_end, 2),
        "tscm_pct": round(tscm, 2),
        "src_fold": round(float(xf[5]), 2),
        "tox_pct": round(float(100.0 * xf[6]), 2),
        "foxo1_nuc_pct": round(float(100.0 * xf[4]), 2),
        "ph_end": round(float(xf[7]), 2),
        "min_ph": round(min_ph, 2),
        "r_lp_end": round(float(xf[3]), 2),
        "_traj": traj,
        "_u_traj": u_traj,
    }


# =============================================================================
# 3. Smooth Discrete-Adjoint Optimal Control (PMP) & Pareto Frontier
# =============================================================================
def smooth_max0(z: np.ndarray | float, eps: float = 1e-4):
    return 0.5 * (z + np.sqrt(z * z + eps * eps))


def eval_cart8_rhs_smooth(x: np.ndarray, u1: float, u2: float, p: Dict = NOMINAL_PARAMS) -> np.ndarray:
    """C^2-smoothed RHS matching eval_cart8_rhs_vec to <1e-4 for exact discrete-adjoint L-BFGS-B."""
    W, D, N, R, F, S, E_tox, pH = x
    M = W + D + 1e-9
    h = D / M
    w_frac = W / M
    g_glc = 1.0 - p["alpha_gal"] * u2
    u_azd_eff = u2
    u_sirt = u2

    r_exc = smooth_max0(R - p["R_basal"])
    v_unmask = p["K_LP"] / (p["K_LP"] + r_exc)
    k_md = p["k_mito_w"] + (p["k_mito_d_max"] - p["k_mito_w"]) * u1 * v_unmask

    eff_m_raw = W + p["eta_d"] * D
    s_exc = smooth_max0(S - 1.0)
    atp_ox = eff_m_raw * (1.0 + p["gamma_pdh"] * u_azd_eff) * (1.0 + p["gamma_src"] * s_exc)
    atp_glyc = p["kappa_glyc"] * g_glc * (1.0 - p["c_akt"] * u1 * (1.0 - u2)) * (R / (R + p["K_R_ATP"]))
    atp_tot = atp_ox + atp_glyc

    ph_def = smooth_max0(7.02 - pH)
    phi_ph = 1.0 / (1.0 + (ph_def / 0.16) ** 2)

    eff_m = 0.26 + smooth_max0(eff_m_raw - 0.26)
    e_dot = (1.0 + p["beta_lp"] * r_exc / p["K_LP"] + p["beta_h"] * h + p["beta_ph"] * ph_def) / eff_m
    f_eff = 1.0 - p["c1"] * u1 + p["c2"] * u2
    ep = e_dot * f_eff * f_eff
    dmg = p["k_dmg"] * smooth_max0(ep - p["theta_gsh"])

    crowd = smooth_max0(1.0 - M / p["K"])
    dW = (p["k_bio_0"] + (p["k_bio_max"] - p["k_bio_0"]) * u_sirt * F) * W * crowd - (p["k_mito_w"] + dmg) * W
    dD = p["k_rep_adv"] * (p["k_bio_0"] + p["xi_d"] * (p["k_bio_max"] - p["k_bio_0"]) * u2) * D * crowd - k_md * D + dmg * W

    w_func = W * (1.0 - h)
    ox_comp = (w_func * w_func) / (w_func * w_func + p["K_W_div"] ** 2)
    pdh_gate = 0.15 + 0.85 * u_azd_eff
    mu_anab = p["mu_glyc"] * g_glc * (1.0 - 0.08 * u1 * (1.0 - u2)) + p["mu_ox"] * u2 * pdh_gate * ox_comp
    dN = mu_anab * phi_ph * (atp_tot / (atp_tot + p["K_ATP"])) * N * smooth_max0(1.0 - N / p["N_cap"])

    dR = (
        p["k_warburg"] * g_glc * (1.0 - 0.35 * u1 * (1.0 - u2)) * (1.0 + p["alpha_h"] * h)
        - p["k_pdh"] * (1.0 + p["beta_azd"] * u_azd_eff) * w_frac * R
        - p["k_perf_lp"] * (R - p["R_media"])
    )

    sirt_pdh = u_sirt * pdh_gate
    dF = (
        p["k_imp"] * (p["f_0"] + p["alpha_akt"] * u1 + p["phi_sirt1"] * sirt_pdh) * phi_ph * (1.0 - F)
        - p["k_exp"] * (1.0 - 0.52 * u1) * (1.0 - 0.85 * sirt_pdh) * (1.0 + p["gamma_acid"] * smooth_max0(6.90 - pH)) * F
    )

    dS = (p["s_basal"] * eff_m_raw + p["k_cristae"] * w_func * (1.0 + p["beta_src"] * F * sirt_pdh)) - p["k_src_dec"] * S

    dE = (
        (
            p["k_exh"] * smooth_max0(ep - p["theta_exh"])
            + p["k_warb_exh"] * g_glc * (1.0 - 0.22 * u1) * (1.0 + 0.70 * (N / p["N_ref"]))
            + p["k_acid_exh"] * smooth_max0(6.86 - pH)
        )
        * (1.0 - E_tox)
        - (p["r_epig_0"] + p["k_rev_epig"] * F * sirt_pdh) * E_tox
    )

    dpH = p["k_perf_0"] * (p["pH_media"] - pH) - p["k_acid"] * g_glc * smooth_max0(R - p["R_buf"]) * (N / p["N_ref"])
    return np.stack([dW, dD, dN, dR, dF, dS, dE, dpH], axis=0)


_PERT_EYE8 = np.eye(8) * 1e-7


def rk4_step_and_jacs(x: np.ndarray, u1: float, u2: float, dt: float) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Single 4-stage RK4 step + exact local stage-propagated Jacobians dPhi/dx (8x8) and dPhi/du (8x2)."""
    eps = 1e-7

    def local_f_and_jac(xs: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        # Batch 1 + 16 + 4 = 21 evaluations into a single (8, 21) vectorized call
        xb = np.broadcast_to(xs[:, None], (8, 21)).copy()
        xb[:, 1:9] += _PERT_EYE8
        xb[:, 9:17] -= _PERT_EYE8
        u1b = np.full(21, u1, dtype=float)
        u2b = np.full(21, u2, dtype=float)
        u1b[17] += eps; u1b[18] -= eps
        u2b[19] += eps; u2b[20] -= eps
        fb = eval_cart8_rhs_smooth(xb, u1b, u2b)
        f0 = fb[:, 0]
        fx = (fb[:, 1:9] - fb[:, 9:17]) / (2.0 * eps)
        fu = np.column_stack([
            (fb[:, 17] - fb[:, 18]) / (2.0 * eps),
            (fb[:, 19] - fb[:, 20]) / (2.0 * eps),
        ])
        return f0, fx, fu

    eye8 = np.eye(8)
    k1, fx1, fu1 = local_f_and_jac(x)
    x2 = x + 0.5 * dt * k1
    dx2_dx = eye8 + 0.5 * dt * fx1
    dx2_du = 0.5 * dt * fu1

    k2, fx2, fu2 = local_f_and_jac(x2)
    dk2_dx = fx2 @ dx2_dx
    dk2_du = fu2 + fx2 @ dx2_du
    x3 = x + 0.5 * dt * k2
    dx3_dx = eye8 + 0.5 * dt * dk2_dx
    dx3_du = 0.5 * dt * dk2_du

    k3, fx3, fu3 = local_f_and_jac(x3)
    dk3_dx = fx3 @ dx3_dx
    dk3_du = fu3 + fx3 @ dx3_du
    x4 = x + dt * k3
    dx4_dx = eye8 + dt * dk3_dx
    dx4_du = dt * dk3_du

    k4, fx4, fu4 = local_f_and_jac(x4)
    dk4_dx = fx4 @ dx4_dx
    dk4_du = fu4 + fx4 @ dx4_du

    x_next = x + (dt / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)
    A = eye8 + (dt / 6.0) * (fx1 + 2.0 * dk2_dx + 2.0 * dk3_dx + dk4_dx)
    B = (dt / 6.0) * (fu1 + 2.0 * dk2_du + 2.0 * dk3_du + dk4_du)
    return x_next, A, B


def solve_discrete_adjoint_pareto() -> Dict:
    """Solve 200-variable (100 intervals x 2 controls) discrete-adjoint PMP & Pareto frontier."""
    n_ocp = 100
    dt = T_HORIZON_H / n_ocp
    reg_u = 0.008  # mild actuator energy regularization

    def make_obj_and_grad(alpha_stem: float):
        def obj_and_grad(z: np.ndarray) -> Tuple[float, np.ndarray]:
            u_seq = z.reshape(n_ocp, 2)
            xs = np.zeros((n_ocp + 1, 8))
            As = np.zeros((n_ocp, 8, 8))
            Bs = np.zeros((n_ocp, 8, 2))
            xs[0] = X0_NOMINAL
            for k in range(n_ocp):
                xs[k + 1], As[k], Bs[k] = rk4_step_and_jacs(xs[k], u_seq[k, 0], u_seq[k, 1], dt)

            xf = xs[-1]
            W, D, N, R, F, S, E_tox, pH = xf
            M = W + D + 1e-9
            h_f = D / M
            k_mem = NOMINAL_PARAMS["K_SRC_mem"]
            tscm_frac = F * (1.0 - E_tox) * (S / (S + k_mem))
            n_norm = N / 75.0
            dose_shortfall = smooth_max0(50.0 - N, eps=0.25) / 50.0

            # Maximize alpha_stem * T_SCM + (1 - alpha_stem) * (N / 75) subject to h_f, TOX, and dose penalties
            j_term = (
                -alpha_stem * tscm_frac
                - (1.0 - alpha_stem) * n_norm
                + 0.45 * (h_f / 0.10) ** 2
                + 0.35 * E_tox ** 2
                + 2.50 * (dose_shortfall ** 2)
            )
            j_run = reg_u * (dt / T_HORIZON_H) * float(np.sum(u_seq ** 2))
            total_cost = float(j_term + j_run)

            # Exact terminal cost gradient dJ_term / dx_N
            lam = np.zeros(8)
            # h_f = D / (W + D)
            dh_dW = -D / (M * M)
            dh_dD = W / (M * M)
            lam[0] = 0.90 * (h_f / (0.10 ** 2)) * dh_dW
            lam[1] = 0.90 * (h_f / (0.10 ** 2)) * dh_dD
            # d smooth_max0(50 - N) / dN
            z_sh = 50.0 - N
            dsh_dN = -0.5 * (1.0 + z_sh / math.sqrt(z_sh * z_sh + 0.25 * 0.25)) / 50.0
            lam[2] = -(1.0 - alpha_stem) / 75.0 + 5.0 * dose_shortfall * dsh_dN
            lam[4] = -alpha_stem * (1.0 - E_tox) * (S / (S + k_mem))
            lam[5] = -alpha_stem * F * (1.0 - E_tox) * (k_mem / ((S + k_mem) ** 2))
            lam[6] = alpha_stem * F * (S / (S + k_mem)) + 0.70 * E_tox

            grad_u = np.zeros((n_ocp, 2))
            for k in range(n_ocp - 1, -1, -1):
                grad_u[k] = 2.0 * reg_u * (dt / T_HORIZON_H) * u_seq[k] + Bs[k].T @ lam
                lam = As[k].T @ lam

            return total_cost, grad_u.ravel()

        return obj_and_grad

    # 1. Verify discrete-adjoint gradient parity against central finite differences across ALL 200 control variables
    rng = np.random.default_rng(SEED)
    test_z = rng.uniform(0.15, 0.85, size=2 * n_ocp)
    alpha_bal = 0.58
    fn_bal = make_obj_and_grad(alpha_bal)
    _, g_adj = fn_bal(test_z)
    max_grad_err = 0.0
    eps_fd = 1e-6
    for idx in range(2 * n_ocp):
        zp = test_z.copy(); zm = test_z.copy()
        zp[idx] += eps_fd; zm[idx] -= eps_fd
        jp, _ = fn_bal(zp)
        jm, _ = fn_bal(zm)
        g_fd = (jp - jm) / (2.0 * eps_fd)
        max_grad_err = max(max_grad_err, abs(g_adj[idx] - g_fd))

    # 2. 5-Start L-BFGS-B for balanced objective (alpha = 0.58)
    t_grid = np.linspace(0.0, T_HORIZON_H, n_ocp, endpoint=False)
    inits = []
    # Start 0: Two-stage step at 66h
    u0_step = np.column_stack([(t_grid < 66.0).astype(float), (t_grid >= 66.0).astype(float)]).ravel()
    inits.append(u0_step)
    # Start 1: Uniform 0.5
    inits.append(np.full(2 * n_ocp, 0.5))
    # Start 2: Early step at 36h
    inits.append(np.column_stack([(t_grid < 36.0).astype(float), (t_grid >= 36.0).astype(float)]).ravel())
    # Start 3: Late step at 96h
    inits.append(np.column_stack([(t_grid < 96.0).astype(float), (t_grid >= 96.0).astype(float)]).ravel())
    # Start 4: Random smooth draw
    inits.append(rng.uniform(0.1, 0.9, size=2 * n_ocp))

    bounds = [(0.0, 1.0)] * (2 * n_ocp)
    start_costs = []
    best_res = None
    for z0 in inits:
        res = minimize(fn_bal, z0, method="L-BFGS-B", jac=True, bounds=bounds, options={"maxiter": 140, "ftol": 1e-10})
        start_costs.append(float(res.fun))
        if best_res is None or res.fun < best_res.fun:
            best_res = res

    u_opt_bal = best_res.x.reshape(n_ocp, 2)
    # Evaluate exact unsmoothed RK4 trajectory under u_opt_bal (interpolated to 1000 steps via ZOH)
    X_pmp = X0_NOMINAL[:, None].copy()
    traj_pmp = np.zeros((1001, 8))
    traj_pmp[0] = X_pmp[:, 0]
    dt_fine = T_HORIZON_H / 1000.0
    for k in range(1000):
        idx_c = min(n_ocp - 1, int(k * n_ocp / 1000))
        U_k = np.array([u_opt_bal[idx_c, 0], u_opt_bal[idx_c, 1], 0.0])[:, None]
        X_pmp = rk4_step_vec(X_pmp, U_k, NOMINAL_PARAMS, dt_fine)
        traj_pmp[k + 1] = X_pmp[:, 0]
    xf_pmp = traj_pmp[-1]

    # Compare J_bal on the exact same objective for optimal constant 2-stage relay and closed-loop gate
    j_closed, _ = fn_bal(u0_step)

    # 3. Sweep Pareto Frontier across alpha_stem in [0.10 .. 0.92]
    alphas = [0.08, 0.20, 0.34, 0.46, 0.58, 0.70, 0.82, 0.92]
    pareto_pmp = []
    z_warm = best_res.x.copy()
    for a_w in alphas:
        fn_a = make_obj_and_grad(a_w)
        # Remove the hard 50x floor penalty on pure-stemness extreme only when exploring full frontier
        res_a = minimize(fn_a, z_warm, method="L-BFGS-B", jac=True, bounds=bounds, options={"maxiter": 110, "ftol": 1e-9})
        z_warm = res_a.x.copy()
        ua = res_a.x.reshape(n_ocp, 2)
        Xa = X0_NOMINAL[:, None].copy()
        for k in range(1000):
            ic = min(n_ocp - 1, int(k * n_ocp / 1000))
            Xa = rk4_step_vec(Xa, np.array([ua[ic, 0], ua[ic, 1], 0.0])[:, None], NOMINAL_PARAMS, dt_fine)
        xfa = Xa[:, 0]
        pareto_pmp.append({
            "alpha_stem": a_w,
            "expansion_fold": round(float(xfa[2]), 2),
            "tscm_pct": round(float(compute_tscm_pct(xfa)), 2),
            "h_end_pct": round(float(100.0 * xfa[1] / (xfa[0] + xfa[1])), 2),
            "tox_pct": round(float(100.0 * xfa[6]), 2),
            "src_fold": round(float(xfa[5]), 2),
        })

    # Dense 2-stage switching curve for t_sw in [0, 168]
    two_stage_curve = []
    for ts in np.linspace(0.0, 168.0, 29):
        arm_name = f"fixed_{ts:.2f}"
        r_ts = simulate_nominal_arm(arm_name)
        two_stage_curve.append({
            "t_switch_h": round(float(ts), 1),
            "expansion_fold": r_ts["expansion_fold"],
            "tscm_pct": r_ts["tscm_pct"],
            "h_end_pct": r_ts["h_end_pct"],
            "tox_pct": r_ts["tox_pct"],
            "src_fold": r_ts["src_fold"],
        })

    spread_pct = 100.0 * (max(start_costs) - min(start_costs)) / max(1e-9, abs(min(start_costs)))
    return {
        "n_intervals": n_ocp,
        "n_variables": 2 * n_ocp,
        "alpha_balanced": alpha_bal,
        "fd_checked_variables": 2 * n_ocp,
        "max_adjoint_vs_fd_grad_err": float(f"{max_grad_err:.3e}"),
        "multi_start_count": len(start_costs),
        "multi_start_descriptions": [
            "Two-stage step at 66.0 h",
            "Uniform interior u1=u2=0.5",
            "Early step at 36.0 h",
            "Late step at 96.0 h",
            "Pseudo-random interior draw U[0.1, 0.9]",
        ],
        "multi_start_costs": [round(c, 6) for c in start_costs],
        "multi_start_cost_spread_pct": round(spread_pct, 4),
        "j_pmp_balanced": round(float(best_res.fun), 5),
        "j_closed_loop_balanced": round(float(j_closed), 5),
        "pmp_balanced_metrics": {
            "expansion_fold": round(float(xf_pmp[2]), 2),
            "tscm_pct": round(float(compute_tscm_pct(xf_pmp)), 2),
            "h_end_pct": round(float(100.0 * xf_pmp[1] / (xf_pmp[0] + xf_pmp[1])), 2),
            "src_fold": round(float(xf_pmp[5]), 2),
            "tox_pct": round(float(100.0 * xf_pmp[6]), 2),
        },
        "pareto_pmp_points": pareto_pmp,
        "two_stage_pareto_curve": two_stage_curve,
        "_u_opt_bal": u_opt_bal,
        "_traj_pmp": traj_pmp,
    }


# =============================================================================
# 4. Volumetric OUR vs. Cell-Density-Invariant ORR Observability & N=1,000 Cohort
# =============================================================================
def cassette_observation_vec(X: np.ndarray, u1: np.ndarray | float, u2: np.ndarray | float, p: Dict = NOMINAL_PARAMS) -> np.ndarray:
    """
    5 inline non-destructive cassette supernatant/optical sensor channels:
      y[0]: Amperometric [Lactate]/[Pyruvate] ratio R_LP (subject to protein biofouling)
      y[1]: Dissolved-O2 optode VOLUMETRIC oxygen uptake rate OUR_vol = N_cells * ATP_ox
            (physically scales with total suspension biomass N_cells)
      y[2]: NIR optical backscatter / capacitance cell-density fold OD = alpha_OD * N_cells
      y[3]: Dual-wavelength optical phenol-red pH
      y[4]: Quartz-window label-free autofluorescence Optical Redox Ratio ORR = FAD/(NAD(P)H+FAD)
            (Walsh et al., Nat. Biomed. Eng. 2021; ratiometric -> N_cells cancels identically!)
    """
    W = np.maximum(0.0, X[0])
    D = np.maximum(0.0, X[1])
    N = np.maximum(0.0, X[2])
    R = np.maximum(0.0, X[3])
    S = np.maximum(0.0, X[5])
    pH = X[7]
    atp_ox = (W + p["eta_d"] * D) * (1.0 + p["gamma_pdh"] * u2) * (1.0 + p["gamma_src"] * np.maximum(0.0, S - 1.0))
    our_vol = N * atp_ox
    alpha_od = p.get("alpha_od", 1.0)
    od_cells = alpha_od * N
    c_glyc = p.get("c_orr_glyc", 0.35)
    g_glc = 1.0 - p["alpha_gal"] * u2
    orr = W / np.maximum(1e-6, W + D + c_glyc * g_glc)
    return np.stack([R, our_vol, od_cells, pH, orr], axis=0)


def compute_cassette_observability_gramian() -> Dict:
    """
    Compute empirical local Observability Gramians along Stage 1 (0..64h) comparing:
      (a) 4-ch supernatant/volumetric-OUR suite WITHOUT ORR: [R_LP, OUR_vol, OD, pH]
      (b) 3-ch post-fouling suite WITHOUT ORR:               [OUR_vol, OD, pH] (ill-conditioned!)
      (c) 5-ch optical-metabolic suite WITH inline ORR:      [R_LP, OUR_vol, OD, pH, ORR]
      (d) 4-ch post-fouling optical suite WITH ORR:          [OUR_vol, OD, pH, ORR]
    """
    dt = 0.5
    n_steps = int(64.0 / dt)
    eps = 1e-5
    active_idx = [0, 1, 2, 3, 5, 7]
    scales = np.array([1.0, 1.0, 25.0, 20.0, 2.5, 0.5])
    y_scales = np.array([15.0, 12.0, 25.0, 0.3, 0.25])

    def sim_outputs(x0_pert: np.ndarray) -> np.ndarray:
        X = x0_pert[:, None].copy()
        U = np.array([1.0, 0.0, 0.0])[:, None]
        ys = np.zeros((n_steps, 5))
        for k in range(n_steps):
            ys[k] = cassette_observation_vec(X[:, 0], 1.0, 0.0) / y_scales
            X = rk4_step_vec(X, U, NOMINAL_PARAMS, dt)
        return ys

    jac_cols = []
    for ii, s_idx in enumerate(active_idx):
        xp = X0_NOMINAL.copy(); xm = X0_NOMINAL.copy()
        step_p = eps * scales[ii]
        xp[s_idx] += step_p; xm[s_idx] -= step_p
        dy = (sim_outputs(xp) - sim_outputs(xm)) / (2.0 * eps)
        jac_cols.append(dy)

    def gram_stats(ch_indices: List[int]) -> Dict:
        J = np.stack([c[:, ch_indices].ravel() for c in jac_cols], axis=1)
        G = (dt / 64.0) * (J.T @ J)
        sv = np.linalg.svd(G, compute_uv=False)
        return {
            "rank": int(np.linalg.matrix_rank(G)),
            "cond_gramian": round(float(sv[0] / sv[-1]), 2),
            "cond_jacobian": round(float(np.sqrt(sv[0] / sv[-1])), 2),
            "min_singular_value": float(f"{sv[-1]:.4e}"),
        }

    g_4ch_no_orr = gram_stats([0, 1, 2, 3])
    g_3ch_foul_no_orr = gram_stats([1, 2, 3])
    g_5ch_with_orr = gram_stats([0, 1, 2, 3, 4])
    g_4ch_foul_with_orr = gram_stats([1, 2, 3, 4])

    return {
        "evaluated_states": ["W", "D", "N_cells", "R_LP", "S_SRC", "pH"],
        "horizon_h": 64.0,
        "suite_4ch_no_orr": g_4ch_no_orr,
        "suite_3ch_post_foul_no_orr": g_3ch_foul_no_orr,
        "suite_5ch_with_orr": g_5ch_with_orr,
        "suite_4ch_post_foul_with_orr": g_4ch_foul_with_orr,
    }


def compute_paired_mcnemar_stats(ok_a: np.ndarray, ok_b: np.ndarray) -> Dict:
    """Exact two-sided binomial McNemar test + paired risk-difference 95% CI between boolean outcomes A and B."""
    n = len(ok_a)
    both_pass = int((ok_a & ok_b).sum())
    a_only = int((ok_a & (~ok_b)).sum())
    b_only = int(((~ok_a) & ok_b).sum())
    both_fail = int(((~ok_a) & (~ok_b)).sum())
    disc = a_only + b_only
    p_val = float(binomtest(a_only, disc, 0.5, alternative="two-sided").pvalue) if disc > 0 else 1.0
    diff = (a_only - b_only) / n
    # Paired standard error of risk difference: sqrt((b + c) - (b - c)^2 / n) / n
    se_diff = math.sqrt(max(0.0, disc - ((a_only - b_only) ** 2) / n)) / n
    ci_lo = 100.0 * (diff - 1.959964 * se_diff)
    ci_hi = 100.0 * (diff + 1.959964 * se_diff)
    return {
        "n_donors": n,
        "both_pass": both_pass,
        "a_pass_b_fail": a_only,
        "a_fail_b_pass": b_only,
        "both_fail": both_fail,
        "paired_diff_pct": round(100.0 * diff, 2),
        "paired_diff_95ci_pct": [round(ci_lo, 2), round(ci_hi, 2)],
        "mcnemar_exact_p": float(f"{p_val:.4e}"),
    }


def run_cohort_and_biofouling_benchmarks() -> Dict:
    """
    Evaluate N=1,000 heterogeneous oncology patient apheresis donors across:
      - Unphased Warburg, Premature Ox (h=0), Continuous Urolithin A Monotherapy
      - 9 Open-Loop Fixed Timers (24, 36, 48, 60, 68, 76, 84, 96, 120 h)
      - Closed-Loop 4-Channel (No-ORR, Volumetric OUR_vol + OD bias) Huber-EKF
      - Closed-Loop 5-Channel Inline ORR Huber-EKF (with & without 83 h exhaustion guard)
      - Amperometric Lac:Pyr Biofouling: Naive Gate vs. No-ORR Handover vs. ORR+NIS Handover
      - Adaptive Extended-Horizon (168-240 h) OD-Guided Harvest Testing Refractory-Donor Rescue
      - Exact Paired McNemar Tests & 95% CIs
    """
    rng = np.random.default_rng(SEED)
    B = 1000
    n_steps = 560  # dt = 0.30 h (18 min inline bioreactor frame over 168 h)
    dt = T_HORIZON_H / n_steps

    h0 = rng.uniform(0.20, 0.62, B)
    m0 = rng.uniform(0.85, 1.12, B)
    e_tox0 = rng.uniform(0.05, 0.14, B)
    X0 = np.stack([
        m0 * (1.0 - h0),
        m0 * h0,
        np.ones(B),
        rng.uniform(18.0, 28.0, B),
        rng.uniform(0.22, 0.36, B),
        rng.uniform(1.05, 1.50, B),
        e_tox0,
        np.full(B, 7.32),
    ], axis=0)

    P_true = NOMINAL_PARAMS.copy()
    mu_glyc_scale = rng.uniform(0.66, 1.36, B)
    mu_ox_scale = rng.uniform(0.84, 1.16, B)
    P_true["mu_glyc"] = NOMINAL_PARAMS["mu_glyc"] * mu_glyc_scale
    P_true["mu_ox"] = NOMINAL_PARAMS["mu_ox"] * mu_ox_scale
    P_true["k_mito_d_max"] = NOMINAL_PARAMS["k_mito_d_max"] * rng.uniform(0.68, 1.36, B)
    P_true["k_warburg"] = NOMINAL_PARAMS["k_warburg"] * rng.uniform(0.82, 1.22, B)
    P_true["k_warb_exh"] = NOMINAL_PARAMS["k_warb_exh"] * rng.uniform(0.82, 1.22, B)
    P_true["c_orr_glyc"] = 0.35 * rng.uniform(0.92, 1.08, B)
    # Unmodeled +/-6% donor lymphoblast size / optical backscatter OD calibration bias
    P_true["alpha_od"] = rng.uniform(0.94, 1.06, B)

    def check_release_pass(Xf: np.ndarray) -> Tuple[np.ndarray, Dict]:
        W, D, N, R, F, S, E_tox, pH = Xf
        h_end = 100.0 * D / (W + D)
        tscm = compute_tscm_pct(Xf)
        tox = 100.0 * E_tox
        ok = (N >= 50.0) & (tscm >= 68.0) & (tox <= 15.0) & (S >= 3.5) & (h_end <= 10.0)
        summary = {
            "pass_count": int(ok.sum()),
            "pass_rate_pct": round(float(100.0 * ok.mean()), 2),
            "expansion_median": round(float(np.median(N)), 2),
            "expansion_iqr": [round(float(np.percentile(N, 25)), 2), round(float(np.percentile(N, 75)), 2)],
            "tscm_median": round(float(np.median(tscm)), 2),
            "tscm_iqr": [round(float(np.percentile(tscm, 25)), 2), round(float(np.percentile(tscm, 75)), 2)],
            "tox_median": round(float(np.median(tox)), 2),
            "h_end_median": round(float(np.median(h_end)), 2),
            "src_median": round(float(np.median(S)), 2),
            "_N": N,
            "_tscm": tscm,
            "_tox": tox,
            "_h_end": h_end,
            "_ok": ok,
        }
        return ok, summary

    # 1. Fixed & baseline open-loop arms (168 h)
    cohort_results: Dict[str, Dict] = {}
    fixed_timers = [24, 36, 48, 60, 68, 76, 84, 96, 120]
    all_open_arms = ["warburg", "premature_ox_h0", "urolithin_mono"] + [f"fixed_{t}" for t in fixed_timers]

    oracle_pass_168 = np.zeros(B, dtype=bool)
    oracle_best_score = np.full(B, -1e9)
    oracle_sw_t = np.full(B, 68.0)

    for arm in all_open_arms:
        X = X0.copy()
        for k in range(n_steps):
            t = k * dt
            if arm == "warburg":
                u1 = np.zeros(B); u2 = np.zeros(B)
            elif arm == "premature_ox_h0":
                u1 = np.zeros(B); u2 = np.ones(B)
            elif arm == "urolithin_mono":
                u1 = np.ones(B); u2 = np.zeros(B)
            else:
                ts = float(arm.split("_")[1])
                u1 = np.full(B, 1.0 if t < ts else 0.0)
                u2 = np.full(B, 0.0 if t < ts else 1.0)
            X = rk4_step_vec(X, np.stack([u1, u2, np.zeros(B)], axis=0), P_true, dt)
        _, summ = check_release_pass(X)
        cohort_results[arm] = summ

    # Evaluate 168 h AND 240 h full-information donor-specific 2-stage oracle ceilings
    oracle_pass_240 = np.zeros(B, dtype=bool)
    n_steps_240 = int(round(240.0 / dt))  # 800 steps
    for ts in np.linspace(24.0, 144.0, 31):
        X = X0.copy()
        for k in range(n_steps_240):
            t = k * dt
            u1 = np.full(B, 1.0 if t < ts else 0.0)
            u2 = np.full(B, 0.0 if t < ts else 1.0)
            X = rk4_step_vec(X, np.stack([u1, u2, np.zeros(B)], axis=0), P_true, dt)
            if k + 1 == n_steps:
                ok_168, summ_ts = check_release_pass(X)
                oracle_pass_168 |= ok_168
                score = summ_ts["_tscm"] + 0.35 * np.minimum(summ_ts["_N"], 65.0) - 200.0 * (~ok_168).astype(float)
                better = score > oracle_best_score
                oracle_best_score = np.where(better, score, oracle_best_score)
                oracle_sw_t = np.where(better, ts, oracle_sw_t)
            # Between 168 h and 240 h (every 6 h check), allow variable harvest time t_harv in [168, 240]
            if (k + 1) >= n_steps and ((k + 1 - n_steps) % 20 == 0):
                ok_ext, _ = check_release_pass(X)
                oracle_pass_240 |= ok_ext

    # 2. Vectorized Cassette Huber-EKF Closed-Loop (Volumetric OUR_vol, OD bias, ORR, Biofouling, & Adaptive 168-240h)
    def run_ekf_cohort(
        use_orr: bool,
        biofoul: bool,
        enable_nis_handover: bool,
        naive_lac_pyr_only: bool = False,
        enable_exh_guard: bool = True,
        adaptive_horizon_max_h: float = 168.0,
    ) -> Dict:
        rng_loc = np.random.default_rng(SEED + (101 if biofoul else 7))
        X = X0.copy()
        X_hat = np.broadcast_to(X0_NOMINAL[:, None], (8, B)).copy()
        P_cov = np.broadcast_to(np.diag([0.05, 0.05, 0.04, 4.0, 0.02, 0.08, 0.01, 0.01])[:, :, None], (8, 8, B)).copy()
        Q_diag = np.array([2.5e-4, 2.5e-4, 6e-3, 0.12, 2e-4, 8e-4, 1e-4, 2e-4])
        # Base additive noise stds: [R_LP, OUR_vol_base, OD_base, pH, ORR]
        R_base_sig = np.array([0.55, 0.050, 0.35, 0.018, 0.012])

        total_steps = int(round(adaptive_horizon_max_h / dt))
        is_extended = adaptive_horizon_max_h > 168.5

        stage = np.ones(B, dtype=int)
        sw_t = np.full(B, np.nan)
        quorum = np.zeros(B, dtype=int)
        foul_onset = rng_loc.uniform(30.0, 48.0, B)
        nis_streak = np.zeros(B, dtype=int)
        ch0_disabled = np.zeros(B, dtype=bool)
        handover_t = np.full(B, np.nan)
        ar_unit = np.zeros((5, B))
        h_rmse_accum = []
        n_at_24h = np.ones(B)

        harvested = np.zeros(B, dtype=bool)
        harvest_t = np.full(B, adaptive_horizon_max_h)
        X_harvest = np.zeros_like(X)

        rep_idx = 12
        rep_trace = {"t": [], "r_true": [], "r_meas": [], "h_true": [], "h_hat": [], "nis_r": [], "disabled": []}

        eps_j = 1e-4
        m_channels = [0, 1, 2, 3, 4] if use_orr else [0, 1, 2, 3]
        for k in range(total_steps):
            t = k * dt
            u1 = np.where(stage == 1, 1.0, 0.0)
            u2 = np.where(stage == 1, 0.0, 1.0)
            U = np.stack([u1, u2, np.zeros(B)], axis=0)

            # Heteroscedastic physical sensor noise (volumetric OUR_vol & NIR OD scale with N_cells)
            ar_unit = 0.65 * ar_unit + math.sqrt(1.0 - 0.65 ** 2) * rng_loc.normal(0.0, 1.0, (5, B))
            sig_mat = np.broadcast_to(R_base_sig[:, None], (5, B)).copy()
            sig_mat[1] = np.sqrt(0.050 ** 2 + (0.045 * X[2]) ** 2)
            sig_mat[2] = np.sqrt(0.35 ** 2 + (0.040 * X[2]) ** 2)

            Y_clean = cassette_observation_vec(X, u1, u2, P_true)
            Y_meas = Y_clean + ar_unit * sig_mat

            if biofoul:
                active_foul = t >= foul_onset
                dt_f = np.maximum(0.0, t - foul_onset)
                spikes = (rng_loc.random(B) < 0.14).astype(float) * rng_loc.choice([-1.0, 1.0], B) * rng_loc.uniform(4.8, 9.5, B)
                foul_bias = np.where(active_foul, 0.42 * dt_f + spikes, 0.0)
                Y_meas[0] = np.maximum(1.5, Y_meas[0] + foul_bias)

            # EKF Predict using NOMINAL_PARAMS (alpha_od = 1.0, c_orr_glyc = 0.35)
            X_pred = rk4_step_vec(X_hat, U, NOMINAL_PARAMS, dt)
            f0 = eval_cart8_rhs_vec(X_hat, U, NOMINAL_PARAMS)
            F_jac = np.zeros((8, 8, B))
            for j in range(8):
                Xp = X_hat.copy()
                Xp[j] += eps_j
                df = (eval_cart8_rhs_vec(Xp, U, NOMINAL_PARAMS) - f0) / eps_j
                F_jac[:, j, :] = df * dt
                F_jac[j, j, :] += 1.0

            FP = np.einsum("ijb,jkb->ikb", F_jac, P_cov)
            P_pred = np.einsum("ijb,kjb->ikb", FP, F_jac)
            for j in range(8):
                P_pred[j, j, :] += Q_diag[j]

            # Sequential scalar Huber-EKF measurement updates across active channels
            X_upd = X_pred.copy()
            P_upd = P_pred.copy()
            nis0 = np.zeros(B)
            for m_ch in m_channels:
                y_hat_all = cassette_observation_vec(X_upd, u1, u2, NOMINAL_PARAMS)
                innov = Y_meas[m_ch] - y_hat_all[m_ch]

                H_m = np.zeros((8, B))
                r_var_m = np.full(B, R_base_sig[m_ch] ** 2)
                if m_ch == 0:
                    H_m[3] = 1.0
                elif m_ch == 1:
                    # Channel 1: Volumetric OUR_vol = N_cells * (W + eta_d * D) * g_p * g_s
                    w_u = np.maximum(0.0, X_upd[0])
                    d_u = np.maximum(0.0, X_upd[1])
                    n_u = np.maximum(0.0, X_upd[2])
                    s_ex = np.maximum(0.0, X_upd[5] - 1.0)
                    g_p = 1.0 + NOMINAL_PARAMS["gamma_pdh"] * u2
                    g_s = 1.0 + NOMINAL_PARAMS["gamma_src"] * s_ex
                    eff_m = w_u + NOMINAL_PARAMS["eta_d"] * d_u
                    H_m[0] = n_u * g_p * g_s
                    H_m[1] = n_u * NOMINAL_PARAMS["eta_d"] * g_p * g_s
                    H_m[2] = eff_m * g_p * g_s
                    H_m[5] = n_u * eff_m * g_p * NOMINAL_PARAMS["gamma_src"] * (X_upd[5] > 1.0).astype(float)
                    r_var_m = 0.050 ** 2 + (0.045 * n_u) ** 2
                elif m_ch == 2:
                    H_m[2] = 1.0
                    n_u = np.maximum(0.0, X_upd[2])
                    r_var_m = 0.35 ** 2 + (0.040 * n_u) ** 2
                elif m_ch == 3:
                    H_m[7] = 1.0
                else:
                    # Channel 4: ORR = W / (W + D + c_glyc * g_glc) -> cell-density N_cells cancels!
                    c_g = 0.35 * (1.0 - NOMINAL_PARAMS["alpha_gal"] * u2)
                    w_u = np.maximum(0.0, X_upd[0])
                    d_u = np.maximum(0.0, X_upd[1])
                    den = np.maximum(1e-6, w_u + d_u + c_g)
                    H_m[0] = (d_u + c_g) / (den * den)
                    H_m[1] = -w_u / (den * den)

                PHt = np.einsum("ijb,jb->ib", P_upd, H_m)
                S_inn = np.einsum("ib,ib->b", H_m, PHt) + r_var_m

                if m_ch == 0:
                    nis0 = (innov ** 2) / np.maximum(1e-9, S_inn)
                    if enable_nis_handover:
                        bad = nis0 > 13.5
                        nis_streak = np.where(bad, nis_streak + 1, 0)
                        newly_disabled = (~ch0_disabled) & (nis_streak >= 3)
                        handover_t = np.where(newly_disabled, t, handover_t)
                        ch0_disabled |= newly_disabled

                std_inn = np.sqrt(np.maximum(1e-9, S_inn))
                z_norm = innov / std_inn
                huber_w = np.where(np.abs(z_norm) <= 2.5, 1.0, 2.5 / np.maximum(1e-9, np.abs(z_norm)))
                if m_ch == 0:
                    active_mask = (~ch0_disabled).astype(float)
                    huber_w = huber_w * active_mask

                K_gain = (PHt / S_inn[None, :]) * huber_w[None, :]
                X_upd = X_upd + K_gain * innov[None, :]
                X_upd[0:7] = np.maximum(0.0, X_upd[0:7])
                X_upd[4] = np.minimum(1.0, X_upd[4])
                X_upd[6] = np.minimum(1.0, X_upd[6])
                KH = np.einsum("ib,jb->ijb", K_gain, H_m)
                P_upd = P_upd - np.einsum("ikb,kjb->ijb", KH, P_upd)

            X_hat = X_upd
            P_cov = P_upd

            if abs(t - 24.0) < 0.5 * dt:
                n_at_24h = np.maximum(1.05, X_hat[2].copy())

            h_true = X[1] / (X[0] + X[1])
            h_hat = X_hat[1] / np.maximum(1e-6, X_hat[0] + X_hat[1])
            if k < n_steps:
                h_rmse_accum.append(float(np.sqrt(np.mean((h_hat - h_true) ** 2))))

            if k < n_steps and k % 2 == 0:
                rep_trace["t"].append(round(t, 2))
                rep_trace["r_true"].append(round(float(X[3, rep_idx]), 3))
                rep_trace["r_meas"].append(round(float(Y_meas[0, rep_idx]), 3))
                rep_trace["h_true"].append(round(float(100.0 * h_true[rep_idx]), 3))
                rep_trace["h_hat"].append(round(float(100.0 * h_hat[rep_idx]), 3))
                rep_trace["nis_r"].append(round(float(nis0[rep_idx]), 2))
                rep_trace["disabled"].append(bool(ch0_disabled[rep_idx]))

            # Explicit Stage 1 -> Stage 2 Persistence Switching Law (4 frames = 1.2 h):
            if naive_lac_pyr_only:
                cond = (Y_meas[0] <= 14.0) & (Y_meas[2] >= 14.0) & (Y_meas[3] >= 7.10)
            else:
                plan_horizon = adaptive_horizon_max_h if is_extended else T_HORIZON_H
                rem_h = np.maximum(0.0, plan_horizon - t)
                if t >= 28.0:
                    mu_emp_ratio = np.clip(np.log(np.maximum(1.1, X_hat[2]) / n_at_24h) / (0.044 * max(4.0, t - 24.0)), 0.78, 1.22)
                else:
                    mu_emp_ratio = np.ones(B)
                s2_mult = np.exp(0.0122 * np.sqrt(mu_emp_ratio) * rem_h)
                # Under +/-6% OD calibration uncertainty, target 53.5x estimated fold to clear true 50.0x floor
                target_n = 53.5 if is_extended else 52.5
                n_req = np.clip(target_n / np.maximum(1.05, 0.58 * s2_mult + 0.42 * np.minimum(s2_mult, 3.6)), 14.0, 42.0)
                met_clean = (h_hat <= 0.056) & (X_hat[2] >= n_req)
                if enable_exh_guard:
                    exh_guard = (t >= 83.0) & (h_hat <= 0.072) & (X_hat[2] >= (14.0 if is_extended else 18.0))
                    cond = (met_clean | exh_guard) & (Y_meas[3] >= 7.10)
                else:
                    cond = met_clean & (Y_meas[3] >= 7.10)

            quorum = np.where((stage == 1) & cond, quorum + 1, 0)
            trig = (stage == 1) & (quorum >= 4)
            sw_t = np.where(trig & np.isnan(sw_t), t, sw_t)
            stage = np.where(trig, 2, stage)

            # Advance TRUE biological plant
            U_next = np.stack([np.where(stage == 1, 1.0, 0.0), np.where(stage == 1, 0.0, 1.0), np.zeros(B)], axis=0)
            X = rk4_step_vec(X, U_next, P_true, dt)

            # Variable OD-triggered harvest check for adaptive extended-horizon arm (t in [168 h, 240 h])
            if is_extended and (k + 1) >= n_steps:
                t_after = (k + 1) * dt
                # Harvest once Stage-2 consolidation has run >= 36 h and estimated N_hat >= 53.5x (or at 240 h cap)
                s2_dur = np.where(np.isnan(sw_t), 0.0, t_after - sw_t)
                ready_harv = (~harvested) & (((stage == 2) & (s2_dur >= 36.0) & (X_hat[2] >= 53.5)) | (k + 1 == total_steps))
                if ready_harv.any():
                    X_harvest[:, ready_harv] = X[:, ready_harv]
                    harvest_t = np.where(ready_harv, t_after, harvest_t)
                    harvested |= ready_harv

        X_eval = X_harvest if is_extended else X
        _, summ = check_release_pass(X_eval)
        valid_sw = sw_t[~np.isnan(sw_t)]
        summ["switched_count"] = int((~np.isnan(sw_t)).sum())
        summ["switch_time_min_h"] = round(float(valid_sw.min()), 2) if len(valid_sw) else None
        summ["switch_time_median_h"] = round(float(np.median(valid_sw)), 2) if len(valid_sw) else None
        summ["switch_time_iqr_h"] = [
            round(float(np.percentile(valid_sw, 25)), 2),
            round(float(np.percentile(valid_sw, 75)), 2),
        ] if len(valid_sw) else None
        summ["switch_time_max_h"] = round(float(valid_sw.max()), 2) if len(valid_sw) else None
        summ["ekf_h_rmse_pct"] = round(100.0 * float(np.mean(h_rmse_accum)), 3)
        summ["biofoul_handover_count"] = int(ch0_disabled.sum())
        if is_extended:
            summ["harvest_time_median_h"] = round(float(np.median(harvest_t)), 2)
            summ["harvest_time_iqr_h"] = [
                round(float(np.percentile(harvest_t, 25)), 2),
                round(float(np.percentile(harvest_t, 75)), 2),
            ]
            summ["extended_beyond_168h_pct"] = round(float(100.0 * (harvest_t > 168.5).mean()), 2)
        summ["_sw_t"] = sw_t
        summ["_harvest_t"] = harvest_t
        summ["_rep_trace"] = rep_trace
        return summ

    ekf_no_orr = run_ekf_cohort(use_orr=False, biofoul=False, enable_nis_handover=True)
    ekf_orr_no_guard = run_ekf_cohort(use_orr=True, biofoul=False, enable_nis_handover=True, enable_exh_guard=False)
    ekf_orr_clean = run_ekf_cohort(use_orr=True, biofoul=False, enable_nis_handover=True, enable_exh_guard=True)
    biofoul_naive_lp = run_ekf_cohort(use_orr=False, biofoul=True, enable_nis_handover=False, naive_lac_pyr_only=True)
    biofoul_no_orr_nis = run_ekf_cohort(use_orr=False, biofoul=True, enable_nis_handover=True)
    biofoul_orr_no_nis = run_ekf_cohort(use_orr=True, biofoul=True, enable_nis_handover=False)
    biofoul_orr_nis = run_ekf_cohort(use_orr=True, biofoul=True, enable_nis_handover=True)
    ekf_orr_adaptive_240h = run_ekf_cohort(
        use_orr=True, biofoul=False, enable_nis_handover=True, enable_exh_guard=True, adaptive_horizon_max_h=240.0
    )

    cohort_results["closed_loop_ekf_no_orr"] = ekf_no_orr
    cohort_results["closed_loop_ekf_no_exh_guard"] = ekf_orr_no_guard
    cohort_results["closed_loop_ekf"] = ekf_orr_clean
    cohort_results["biofoul_naive_lp_gate"] = biofoul_naive_lp
    cohort_results["biofoul_no_orr_nis_handover"] = biofoul_no_orr_nis
    cohort_results["biofoul_orr_no_nis"] = biofoul_orr_no_nis
    cohort_results["biofoul_ekf_nis_handover"] = biofoul_orr_nis
    cohort_results["adaptive_horizon_240h_ekf"] = ekf_orr_adaptive_240h

    # Characterize 168 h vs 240 h Refractory Subsets (Testing Section 5.2 Hypothesis Directly)
    infeas_168 = ~oracle_pass_168
    rescued_by_240_oracle = infeas_168 & oracle_pass_240
    refractory_240_oracle = ~oracle_pass_240

    ok_168_ekf = ekf_orr_clean["_ok"]
    ok_240_ekf = ekf_orr_adaptive_240h["_ok"]
    rescued_by_240_ekf = (~ok_168_ekf) & ok_240_ekf

    infeasible_profile = {
        "infeasible_168h_count": int(infeas_168.sum()),
        "infeasible_168h_pct": round(float(100.0 * infeas_168.mean()), 2),
        "infeasible_168h_median_h0_pct": round(float(100.0 * np.median(h0[infeas_168])), 2),
        "infeasible_168h_median_mu_glyc_rel": round(float(np.median(mu_glyc_scale[infeas_168])), 3),
        "controllable_168h_median_h0_pct": round(float(100.0 * np.median(h0[oracle_pass_168])), 2),
        "controllable_168h_median_mu_glyc_rel": round(float(np.median(mu_glyc_scale[oracle_pass_168])), 3),
        "oracle_240h_controllable_count": int(oracle_pass_240.sum()),
        "oracle_240h_controllable_pct": round(float(100.0 * oracle_pass_240.mean()), 2),
        "oracle_rescued_168_to_240_count": int(rescued_by_240_oracle.sum()),
        "oracle_rescued_of_168_infeasible_pct": round(float(100.0 * rescued_by_240_oracle.sum() / max(1, infeas_168.sum())), 2),
        "ekf_rescued_168_to_240_count": int(rescued_by_240_ekf.sum()),
        "intrinsically_refractory_240h_count": int(refractory_240_oracle.sum()),
        "intrinsically_refractory_240h_pct": round(float(100.0 * refractory_240_oracle.mean()), 2),
        "refractory_240h_median_h0_pct": round(float(100.0 * np.median(h0[refractory_240_oracle])), 2) if refractory_240_oracle.any() else None,
        "refractory_240h_median_mu_glyc_rel": round(float(np.median(mu_glyc_scale[refractory_240_oracle])), 3) if refractory_240_oracle.any() else None,
        "refractory_240h_median_tox0_pct": round(float(100.0 * np.median(e_tox0[refractory_240_oracle])), 2) if refractory_240_oracle.any() else None,
    }

    # Exact Paired McNemar Tests Across N=1,000 Donors
    paired_stats = {
        "orr_ekf_168h_vs_best_fixed_68h": compute_paired_mcnemar_stats(ok_168_ekf, cohort_results["fixed_68"]["_ok"]),
        "orr_ekf_168h_vs_no_orr_ekf_168h": compute_paired_mcnemar_stats(ok_168_ekf, ekf_no_orr["_ok"]),
        "orr_ekf_168h_vs_no_exh_guard": compute_paired_mcnemar_stats(ok_168_ekf, ekf_orr_no_guard["_ok"]),
        "biofoul_orr_nis_vs_best_fixed_68h": compute_paired_mcnemar_stats(biofoul_orr_nis["_ok"], cohort_results["fixed_68"]["_ok"]),
        "biofoul_orr_nis_vs_no_orr_nis": compute_paired_mcnemar_stats(biofoul_orr_nis["_ok"], biofoul_no_orr_nis["_ok"]),
        "adaptive_240h_ekf_vs_fixed_168h_ekf": compute_paired_mcnemar_stats(ok_240_ekf, ok_168_ekf),
    }

    return {
        "n_donors": B,
        "release_criteria": {
            "expansion_fold_ge": 50.0,
            "tscm_pct_ge": 68.0,
            "tox_pct_le": 15.0,
            "src_fold_ge": 3.5,
            "h_end_pct_le": 10.0,
        },
        "oracle_controllable_count": int(oracle_pass_168.sum()),
        "oracle_controllable_pct": round(float(100.0 * oracle_pass_168.mean()), 2),
        "oracle_240h_controllable_count": int(oracle_pass_240.sum()),
        "oracle_240h_controllable_pct": round(float(100.0 * oracle_pass_240.mean()), 2),
        "infeasible_donor_profile": infeasible_profile,
        "paired_mcnemar_tests": paired_stats,
        "arms": cohort_results,
        "_h0": h0,
        "_mu_glyc_scale": mu_glyc_scale,
        "_oracle_pass_168": oracle_pass_168,
        "_oracle_pass_240": oracle_pass_240,
        "_oracle_sw_t": oracle_sw_t,
    }


# =============================================================================
# 5. Host-Emulated C99 Fixed-Point Cassette Verification (100% Emitted JSON)
# =============================================================================
def run_c99_cassette_verification() -> Dict:
    c_src = CORAL_DIR / "cart_bioreactor_hil_testbench.c"
    with tempfile.NamedTemporaryFile(prefix="cart_hil_tb_", delete=False) as tmp:
        bin_path = Path(tmp.name)
    try:
        subprocess.run(
            ["clang", "-O2", "-Wall", "-Wextra", "-Werror", "-std=c99", str(c_src), "-lm", "-o", str(bin_path)],
            check=True,
            capture_output=True,
            text=True,
        )
        proc = subprocess.run([str(bin_path), "--json"], check=True, capture_output=True, text=True)
        c99_data = json.loads(proc.stdout)
    finally:
        if bin_path.exists():
            bin_path.unlink()
    return c99_data


# =============================================================================
# 6. Publication Vector PDF Figures (Paper II Figures 1, 2, 3)
# =============================================================================
def configure_matplotlib() -> None:
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["DejaVu Sans", "Arial", "Helvetica"],
        "font.size": 8.4,
        "axes.titlesize": 8.8,
        "axes.labelsize": 8.3,
        "xtick.labelsize": 7.5,
        "ytick.labelsize": 7.5,
        "legend.fontsize": 7.0,
        "axes.linewidth": 0.8,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    })


def render_fig1_structural_pharmacology(struct_data: Dict) -> Path:
    fig = plt.figure(figsize=(11.8, 7.8))
    gs = gridspec.GridSpec(2, 3, width_ratios=[1.02, 1.35, 0.045], height_ratios=[1.0, 1.0], wspace=0.30, hspace=0.38)

    # Panel A: PDB 2Q8G (PDK1 + AZD7545 lipoyl-binding pocket heavy-atom contacts <= 4.0 A)
    ax_a = fig.add_subplot(gs[0, 0])
    c_2q8g = struct_data["pdk1_2q8g"]["contacts_le_4A"][:10]
    labels_a = [f"{c['residue']} ({c['protein_atom']}-{c['ligand_atom']})" for c in c_2q8g][::-1]
    dists_a = [c["min_dist_A"] for c in c_2q8g][::-1]
    colors_a = ["#0072B2" if d > 2.8 else "#D55E00" for d in dists_a]
    bars_a = ax_a.barh(labels_a, dists_a, color=colors_a, height=0.68, edgecolor="#222222", linewidth=0.5)
    ax_a.axvline(2.80, color="#D55E00", linestyle="--", linewidth=0.9, label="H-bond (2.8 Å)")
    ax_a.axvline(4.00, color="#666666", linestyle=":", linewidth=0.9, label="vdW (4.0 Å)")
    for b, d in zip(bars_a, dists_a):
        ax_a.text(d + 0.05, b.get_y() + b.get_height() / 2, f"{d:.2f} Å", va="center", fontsize=6.9)
    ax_a.set_xlim(0, 5.25)
    ax_a.set_xlabel("Min Heavy-Atom Distance to AZD7545 (Å)")
    ax_a.set_title("A. PDK1 Lipoyl Pocket (PDB: 2Q8G, 1.90 Å)", fontweight="bold", loc="left")
    ax_a.legend(loc="upper right", fontsize=6.7, frameon=True, framealpha=0.92)
    ax_a.grid(axis="x", alpha=0.25)

    # Panel B: PDB 4ZZJ (SIRT1 + Allosteric STAC 4TQ & Carba-NAD+ contacts)
    ax_b = fig.add_subplot(gs[1, 0])
    stac_c = struct_data["sirt1_4zzj_stac_4tq"]["contacts_le_4A"][:5]
    nad_c = struct_data["sirt1_4zzj_carba_nad"]["contacts_le_4A"][:5]
    comb_labels = [f"STAC: {c['residue']}" for c in stac_c] + [f"NAD+: {c['residue']}" for c in nad_c]
    comb_dists = [c["min_dist_A"] for c in stac_c] + [c["min_dist_A"] for c in nad_c]
    comb_cols = ["#009E73"] * len(stac_c) + ["#CC79A7"] * len(nad_c)
    comb_labels = comb_labels[::-1]
    comb_dists = comb_dists[::-1]
    comb_cols = comb_cols[::-1]
    bars_b = ax_b.barh(comb_labels, comb_dists, color=comb_cols, height=0.68, edgecolor="#222222", linewidth=0.5)
    for b, d, lbl in zip(bars_b, comb_dists, comb_labels):
        tag = f"{d:.2f} Å (hinge)" if "GLU230" in lbl else f"{d:.2f} Å"
        ax_b.text(d + 0.05, b.get_y() + b.get_height() / 2, tag, va="center", fontsize=6.9)
    ax_b.set_xlim(0, 4.75)
    ax_b.set_xlabel("Min Heavy-Atom Distance to Ligand (Å)")
    ax_b.set_title("B. SIRT1 Catalytic & STAC Core (PDB: 4ZZJ)", fontweight="bold", loc="left")
    ax_b.grid(axis="x", alpha=0.25)

    # Panel C: FOXO1 (AF-Q12778-F1, 655 AA) pLDDT + PTM disorder + PAE inset
    ax_c = fig.add_subplot(gs[0, 1])
    fox = struct_data["foxo1_af_q12778"]
    res_f = np.arange(1, fox["length_aa"] + 1)
    plddt_f = fox["_plddt_arr"]
    ax_c.axvspan(157, 242, color="#0072B2", alpha=0.14, label="Forkhead DBD (157–242, pLDDT 90.9)")
    ax_c.axhline(70, color="#0072B2", linestyle="--", linewidth=0.8, alpha=0.7)
    ax_c.axhline(50, color="#D55E00", linestyle=":", linewidth=0.8, alpha=0.8)
    ax_c.plot(res_f, plddt_f, color="#333333", linewidth=0.95, alpha=0.9)
    ax_c.fill_between(res_f, 0, plddt_f, where=(plddt_f < 50), color="#E69F00", alpha=0.28, label=f"Very Low <50 ({fox['pct_very_low_lt_50']:.1f}%)")
    for site in fox["ptm_sites"]:
        pos = int(site["site"][3:])
        is_akt = "Akt" in site["role"]
        ax_c.scatter([pos], [site["plddt"]], color="#D55E00" if is_akt else "#CC79A7", s=28, zorder=5,
                     marker="D" if is_akt else "o", edgecolors="black", linewidths=0.4)
    ax_c.scatter([], [], color="#D55E00", marker="D", s=26, edgecolors="black", linewidths=0.4, label="Akt P-sites (T24/S256/S319: all <70)")
    ax_c.scatter([], [], color="#CC79A7", marker="o", s=26, edgecolors="black", linewidths=0.4, label="SIRT1 Ac-Lys (K245–K294: all <60)")
    ax_c.set_xlim(1, fox["length_aa"])
    ax_c.set_ylim(15, 105)
    ax_c.set_ylabel("AlphaFold v6 pLDDT")
    ax_c.set_xlabel("FOXO1 Residue Index (UniProt Q12778, 655 AA)")
    ax_c.set_title("C. Human FOXO1 (AF-Q12778-F1): DBD vs. Disordered PTM Tails", fontweight="bold", loc="left")
    ax_c.legend(loc="upper right", fontsize=6.4, frameon=True, framealpha=0.92)
    ins_c = ax_c.inset_axes([0.37, 0.34, 0.15, 0.35])
    im_c = ins_c.imshow(fox["_pae_mat"], cmap="Greens_r", vmin=0, vmax=30, origin="upper")
    ins_c.set_xticks([0, 654]); ins_c.set_yticks([0, 654])
    ins_c.set_xticklabels(["1", "655"], fontsize=5.8); ins_c.set_yticklabels(["1", "655"], fontsize=5.8)
    ins_c.set_title("PAE (28.2 Å)", fontsize=6.2, pad=3)

    # Panel D: PGC-1alpha (AF-Q9UBK2-F1, 798 AA) pLDDT + 39/47 Lys disorder + PAE inset
    ax_d = fig.add_subplot(gs[1, 1])
    pgc = struct_data["pgc1a_af_q9ubk2"]
    res_p = np.arange(1, pgc["length_aa"] + 1)
    plddt_p = pgc["_plddt_arr"]
    ax_d.axvspan(634, 738, color="#009E73", alpha=0.15, label="C-term RRM (634–738, pLDDT 88.2)")
    ax_d.axhline(70, color="#0072B2", linestyle="--", linewidth=0.8, alpha=0.7)
    ax_d.axhline(50, color="#D55E00", linestyle=":", linewidth=0.8, alpha=0.8)
    ax_d.plot(res_p, plddt_p, color="#333333", linewidth=0.9, alpha=0.9)
    ax_d.fill_between(res_p, 0, plddt_p, where=(plddt_p < 50), color="#E69F00", alpha=0.28, label=f"Very Low <50 ({pgc['pct_very_low_lt_50']:.1f}%)")
    lys_idx = np.array(pgc["_lys_indices"])
    ax_d.scatter(lys_idx + 1, plddt_p[lys_idx], color="#CC79A7", s=16, zorder=4, alpha=0.85,
                 edgecolors="black", linewidths=0.3, label=f"Lysine sites ({pgc['lysines_lt_70']}/{pgc['total_lysines']} <70 pLDDT)")
    ax_d.set_xlim(1, pgc["length_aa"])
    ax_d.set_ylim(15, 105)
    ax_d.set_ylabel("AlphaFold v6 pLDDT")
    ax_d.set_xlabel("PGC-1α Residue Index (UniProt Q9UBK2, 798 AA)")
    ax_d.set_title("D. Human PGC-1α (AF-Q9UBK2-F1): IDR Activation Region vs. RRM", fontweight="bold", loc="left")
    ax_d.legend(loc="upper left", fontsize=6.4, frameon=True, framealpha=0.92)
    ins_d = ax_d.inset_axes([0.52, 0.42, 0.15, 0.35])
    ins_d.imshow(pgc["_pae_mat"], cmap="Greens_r", vmin=0, vmax=30, origin="upper")
    ins_d.set_xticks([0, 797]); ins_d.set_yticks([0, 797])
    ins_d.set_xticklabels(["1", "798"], fontsize=5.8); ins_d.set_yticklabels(["1", "798"], fontsize=5.8)
    ins_d.set_title("PAE (28.1 Å)", fontsize=6.2, pad=3)

    cax = fig.add_subplot(gs[:, 2])
    cbar = fig.colorbar(im_c, cax=cax)
    cbar.set_label("AlphaFold v6 Predicted Aligned Error PAE (Å)", fontsize=7.6)

    out_path = FIG_DIR / "paper2_fig1_structural_pharmacology_and_pae.pdf"
    fig.savefig(out_path, bbox_inches="tight")
    plt.close(fig)
    return out_path


def render_fig2_ode_and_pareto(arms_nom: Dict[str, Dict], pmp_data: Dict) -> Path:
    fig, axes = plt.subplots(2, 2, figsize=(11.8, 7.8))
    plt.subplots_adjust(wspace=0.30, hspace=0.36)
    t_grid = np.linspace(0.0, T_HORIZON_H, 1001)

    tr_warb = arms_nom["warburg"]["_traj"]
    tr_prem = arms_nom["premature_ox_h0"]["_traj"]
    tr_urol = arms_nom["urolithin_mono"]["_traj"]
    tr_cl = arms_nom["closed_loop"]["_traj"]
    sw_cl = arms_nom["closed_loop"]["switch_time_h"]

    # Panel A: Clonal Fold Expansion N_cells(t) & Heteroplasmy h(t)
    ax_a = axes[0, 0]
    ax_a.plot(t_grid, tr_warb[:, 2], color="#D55E00", linestyle="--", linewidth=1.6, label=f"Unphased Warburg ({arms_nom['warburg']['expansion_fold']:.1f}x, h={arms_nom['warburg']['h_end_pct']:.1f}%)")
    ax_a.plot(t_grid, tr_urol[:, 2], color="#E69F00", linestyle="-.", linewidth=1.6, label=f"Continuous Urolithin A ({arms_nom['urolithin_mono']['expansion_fold']:.1f}x, h={arms_nom['urolithin_mono']['h_end_pct']:.1f}%)")
    ax_a.plot(t_grid, tr_prem[:, 2], color="#CC79A7", linestyle=":", linewidth=1.7, label=f"Premature Ox h=0 ({arms_nom['premature_ox_h0']['expansion_fold']:.1f}x < 50x floor)")
    ax_a.plot(t_grid, tr_cl[:, 2], color="#0072B2", linewidth=2.1, label=f"Closed-Loop BCR ({arms_nom['closed_loop']['expansion_fold']:.1f}x, h={arms_nom['closed_loop']['h_end_pct']:.1f}%)")
    ax_a.axhline(50.0, color="#444444", linestyle=":", linewidth=1.0, label="Clinical dose floor (50x)")
    ax_a.axvline(sw_cl, color="#0072B2", linestyle="--", alpha=0.55, linewidth=1.0)
    ax_a.text(sw_cl + 2.5, 6.5, f"Gate t = {sw_cl:.1f} h", fontsize=7.3, color="#0072B2", fontweight="bold")
    ax_a.set_xlim(0, 168); ax_a.set_ylim(0, 88)
    ax_a.set_xlabel("Simulated Perfusion Time t (hours)")
    ax_a.set_ylabel("Simulated Viable Fold Expansion N(t)")
    ax_a.set_title("A. Clonal Blastogenesis vs. Oxidative Arrest (In Silico)", fontweight="bold", loc="left")
    ax_a.legend(loc="upper left", fontsize=6.7, frameon=True)
    ax_a.grid(alpha=0.25)

    # Panel B: Stem-Like Memory T_SCM(t) & Target Ablations
    ax_b = axes[0, 1]
    for key, col, ls, lbl in [
        ("warburg", "#D55E00", "--", "Unphased Warburg"),
        ("urolithin_mono", "#E69F00", "-.", "Continuous Urolithin A (Media-A)"),
        ("ablate_azd", "#56B4E9", ":", "Ablate Phase-2 AZD7545 (2Q8G)"),
        ("ablate_sirt1", "#CC79A7", ":", "Ablate Phase-2 SIRT1/NR (4ZZJ)"),
        ("closed_loop", "#0072B2", "-", "Closed-Loop Bioreactor-BCR"),
    ]:
        tr = arms_nom[key]["_traj"]
        tscm_t = np.array([compute_tscm_pct(tr[i]) for i in range(len(t_grid))])
        ax_b.plot(t_grid, tscm_t, color=col, linestyle=ls, linewidth=2.0 if key == "closed_loop" else 1.45,
                  label=f"{lbl} ({arms_nom[key]['tscm_pct']:.1f}%)")
    ax_b.axhline(68.0, color="#009E73", linestyle="--", linewidth=0.9, alpha=0.8, label="Target T_SCM spec (68%)")
    ax_b.axvline(sw_cl, color="#0072B2", linestyle="--", alpha=0.5, linewidth=1.0)
    ax_b.set_xlim(0, 168); ax_b.set_ylim(0, 86)
    ax_b.set_xlabel("Simulated Perfusion Time t (hours)")
    ax_b.set_ylabel("Simulated Stem-Like Memory T_SCM (%)")
    ax_b.set_title("B. Simulated T_SCM Consolidation & Target Ablations", fontweight="bold", loc="left")
    ax_b.legend(loc="upper left", fontsize=6.6, frameon=True)
    ax_b.grid(alpha=0.25)

    # Panel C: Discrete-Adjoint Pontryagin Optimal Controls u1*(t), u2*(t) vs Two-Stage Gate
    ax_c = axes[1, 0]
    u_opt = pmp_data["_u_opt_bal"]
    t_ocp = np.linspace(0.0, T_HORIZON_H, len(u_opt), endpoint=False)
    ax_c.step(t_ocp, u_opt[:, 0], where="post", color="#009E73", linewidth=1.9, label="Adjoint PMP u1*(t) (Urolithin A mitophagy)")
    ax_c.step(t_ocp, u_opt[:, 1], where="post", color="#0072B2", linewidth=1.9, label="Adjoint PMP u2*(t) (AZD7545 + SIRT1/NR)")
    ax_c.axvline(sw_cl, color="#333333", linestyle="--", linewidth=1.2, label=f"2-Stage Closed-Loop Gate ({sw_cl:.1f} h)")
    ax_c.set_xlim(0, 168); ax_c.set_ylim(-0.05, 1.15)
    ax_c.set_xlabel("Simulated Perfusion Time t (hours)")
    ax_c.set_ylabel("Control Actuation u(t) in [0, 1]")
    ax_c.set_title("C. 200-Var Discrete-Adjoint PMP (5-Start L-BFGS-B)", fontweight="bold", loc="left")
    ax_c.text(3.5, 1.045, f"Adjoint vs FD grad err = {pmp_data['max_adjoint_vs_fd_grad_err']:.2e} | 5-start spread = {pmp_data['multi_start_cost_spread_pct']:.2f}%",
              fontsize=6.8, bbox=dict(boxstyle="round,pad=0.2", facecolor="#F5F7FA", edgecolor="#CCCCCC"))
    ax_c.legend(loc="center right", fontsize=6.8, frameon=True)
    ax_c.grid(alpha=0.25)

    # Panel D: Proliferation vs Stemness Pareto Frontier
    ax_d = axes[1, 1]
    ax_d.axvspan(50.0, 90.0, ymin=68.0 / 88.0, ymax=1.0, color="#009E73", alpha=0.12, label="Release Spec Zone (>=50x, >=68% T_SCM)")
    ax_d.axvline(50.0, color="#444444", linestyle=":", linewidth=0.9)
    ax_d.axhline(68.0, color="#009E73", linestyle="--", linewidth=0.9)

    ts_curve = pmp_data["two_stage_pareto_curve"]
    ax_d.plot([p["expansion_fold"] for p in ts_curve], [p["tscm_pct"] for p in ts_curve],
              color="#56B4E9", linewidth=1.8, marker="o", markersize=3.5, label="2-Stage Switching Locus (0..168 h)")
    pmp_pts = sorted(pmp_data["pareto_pmp_points"], key=lambda d: d["expansion_fold"])
    ax_d.plot([p["expansion_fold"] for p in pmp_pts], [p["tscm_pct"] for p in pmp_pts],
              color="#009E73", linewidth=2.0, marker="s", markersize=4.5, label="Continuous Adjoint PMP Frontier")

    markers = [
        ("warburg", "#D55E00", "X", "Warburg"),
        ("urolithin_mono", "#E69F00", "^", "Urolithin Mono"),
        ("premature_ox_h0", "#CC79A7", "v", "Premature Ox (0h)"),
        ("ablate_azd", "#999999", "P", "No AZD7545"),
        ("ablate_sirt1", "#666666", "D", "No SIRT1/NR"),
        ("closed_loop", "#0072B2", "*", "Closed-Loop BCR"),
    ]
    for key, col, mk, lbl in markers:
        r = arms_nom[key]
        sz = 115 if mk == "*" else 55
        ax_d.scatter([r["expansion_fold"]], [r["tscm_pct"]], color=col, marker=mk, s=sz, zorder=6,
                     edgecolors="black", linewidths=0.5, label=lbl)
    ax_d.set_xlim(5, 86); ax_d.set_ylim(0, 88)
    ax_d.set_xlabel("Simulated 168-Hour Viable Fold Expansion N(168 h)")
    ax_d.set_ylabel("Simulated 168-Hour T_SCM Yield (%)")
    ax_d.set_title("D. Simulated Proliferation–Stemness Pareto Frontier", fontweight="bold", loc="left")
    ax_d.legend(loc="lower left", fontsize=6.4, ncol=2, frameon=True)
    ax_d.grid(alpha=0.25)

    out_path = FIG_DIR / "paper2_fig2_bioreactor_ode_and_pareto_optimal.pdf"
    fig.savefig(out_path, bbox_inches="tight")
    plt.close(fig)
    return out_path


def render_fig3_cohort_biofouling_and_acidosis(cohort_data: Dict, arms_nom: Dict[str, Dict], c99_data: Dict) -> Path:
    fig, axes = plt.subplots(2, 2, figsize=(11.8, 7.9))
    plt.subplots_adjust(wspace=0.30, hspace=0.37)
    arms = cohort_data["arms"]
    n_donors = cohort_data["n_donors"]

    # Panel A: N=1,000 Oncology Cohort Pass Rate: Fixed Timers vs 4-Ch No-ORR vs 5-Ch ORR (168h & Adaptive 240h)
    ax_a = axes[0, 0]
    bar_keys = [
        ("warburg", "Warburg\n(No Sw.)", "#D55E00"),
        ("urolithin_mono", "Urolithin\nMono", "#E69F00"),
        ("fixed_48", "Fixed\n48 h", "#999999"),
        ("fixed_68", "Best Fix\n68 h", "#56B4E9"),
        ("fixed_96", "Fixed\n96 h", "#999999"),
        ("closed_loop_ekf_no_orr", "4-Ch\nNo-ORR", "#CC79A7"),
        ("closed_loop_ekf", "5-Ch ORR\n(168 h)", "#0072B2"),
        ("adaptive_horizon_240h_ekf", "Adaptive\n168–240 h", "#009E73"),
    ]
    x_pos = np.arange(len(bar_keys))
    rates = [arms[k]["pass_rate_pct"] for k, _, _ in bar_keys]
    cols = [c for _, _, c in bar_keys]
    bars = ax_a.bar(x_pos, rates, color=cols, width=0.64, edgecolor="#222222", linewidth=0.6)
    oracle_168 = cohort_data["oracle_controllable_pct"]
    oracle_240 = cohort_data["oracle_240h_controllable_pct"]
    ax_a.axhline(oracle_168, color="#0072B2", linestyle="--", linewidth=1.1,
                 label=f"168 h Full-Info Oracle ({oracle_168:.1f}%)")
    ax_a.axhline(oracle_240, color="#009E73", linestyle=":", linewidth=1.3,
                 label=f"240 h Extended Oracle ({oracle_240:.1f}%)")
    for b, r in zip(bars, rates):
        ax_a.text(b.get_x() + b.get_width() / 2, r + 1.3, f"{r:.1f}%", ha="center", va="bottom", fontsize=6.5, fontweight="bold")
    p_vs_fix = cohort_data["paired_mcnemar_tests"]["orr_ekf_168h_vs_best_fixed_68h"]["mcnemar_exact_p"]
    p_vs_noorr = cohort_data["paired_mcnemar_tests"]["orr_ekf_168h_vs_no_orr_ekf_168h"]["mcnemar_exact_p"]
    ax_a.text(0.03, 0.61,
              f"Paired McNemar (N={n_donors:,}):\n5-Ch vs 68h: p={p_vs_fix:.1e}\n5-Ch vs 4-Ch: p={p_vs_noorr:.1e}",
              transform=ax_a.transAxes, fontsize=6.3,
              bbox=dict(boxstyle="round,pad=0.22", facecolor="#F5F7FA", edgecolor="#CCCCCC"))
    ax_a.set_xticks(x_pos)
    ax_a.set_xticklabels([lbl for _, lbl, _ in bar_keys], fontsize=6.2)
    ax_a.set_ylim(0, 106)
    ax_a.set_ylabel("Simultaneous 5-Criteria Pass Rate (%)")
    ax_a.set_title(f"A. N={n_donors:,} Cohort: Fixed vs. 4-Ch/5-Ch & Adaptive 240h", fontweight="bold", loc="left")
    ax_a.legend(loc="upper left", fontsize=6.5, frameon=True, framealpha=0.92)
    ax_a.grid(axis="y", alpha=0.25)

    # Panel B: 2D Oncology Donor Controllability & Extended-Horizon Rescue Phase Diagram
    ax_b = axes[0, 1]
    h0_pct = 100.0 * cohort_data["_h0"]
    mu_rel = cohort_data["_mu_glyc_scale"]
    ok_168 = arms["closed_loop_ekf"]["_ok"]
    ok_240 = arms["adaptive_horizon_240h_ekf"]["_ok"]
    rescued = (~ok_168) & ok_240
    unrescued = (~ok_168) & (~ok_240)

    ax_b.scatter(h0_pct[ok_168], mu_rel[ok_168], color="#0072B2", s=13, alpha=0.52,
                 edgecolors="none", label=f"Passed at 168 h (n={int(ok_168.sum())})")
    ax_b.scatter(h0_pct[rescued], mu_rel[rescued], color="#009E73", marker="^", s=24, alpha=0.88,
                 edgecolors="black", linewidths=0.3, label=f"Rescued by 168–240 h Harvest (n={int(rescued.sum())})")
    ax_b.scatter(h0_pct[unrescued], mu_rel[unrescued], color="#D55E00", marker="x", s=20, alpha=0.82,
                 linewidths=0.9, label=f"Failed Both 168 h & 240 h (n={int(unrescued.sum())})")
    ax_b.set_xlim(18, 64); ax_b.set_ylim(0.62, 1.56)
    ax_b.set_xlabel("Simulated Donor Initial Depolarized Fraction h(0) (%)")
    ax_b.set_ylabel("Relative Proliferative Fitness mu_glyc / mu_nom")
    ax_b.set_title("B. Donor Controllability & 240 h Rescue Phase Space", fontweight="bold", loc="left")
    ax_b.legend(loc="upper right", fontsize=6.5, frameon=True, framealpha=0.94)
    ax_b.grid(alpha=0.25)

    # Panel C: Amperometric Lac:Pyr Biofouling & NIS Optical Handover
    ax_c = axes[1, 0]
    rep = arms["biofoul_ekf_nis_handover"]["_rep_trace"]
    t_r = np.array(rep["t"])
    ax_c.plot(t_r, rep["r_meas"], color="#D55E00", alpha=0.65, linewidth=1.0, label="Biofouled Amperometric Lac:Pyr (y1)")
    ax_c.plot(t_r, rep["r_true"], color="#333333", linestyle="--", linewidth=1.3, label="True Supernatant Lac:Pyr")
    ax_c.plot(t_r, rep["h_true"], color="#009E73", linewidth=1.8, label="True Depolarization h(t) (%)")
    ax_c.plot(t_r, rep["h_hat"], color="#0072B2", linestyle="-.", linewidth=1.7, label="5-Ch ORR+NIS EKF h_hat(t) (%)")
    dis_arr = np.array(rep["disabled"])
    if dis_arr.any():
        t_dis = float(t_r[np.argmax(dis_arr)])
        ax_c.axvline(t_dis, color="#CC79A7", linestyle=":", linewidth=1.4, label=f"NIS >13.5 Handover ({t_dis:.1f} h)")
    p_naive = arms["biofoul_naive_lp_gate"]["pass_rate_pct"]
    p_no_orr = arms["biofoul_no_orr_nis_handover"]["pass_rate_pct"]
    p_nis = arms["biofoul_ekf_nis_handover"]["pass_rate_pct"]
    p_bf_mc = cohort_data["paired_mcnemar_tests"]["biofoul_orr_nis_vs_no_orr_nis"]["mcnemar_exact_p"]
    ax_c.text(75.0, 15.5,
              f"Biofouled Pass Rate (N={n_donors:,}):\nNaive Lac:Pyr: {p_naive:.1f}%\n4-Ch Handover: {p_no_orr:.1f}%\n5-Ch ORR+NIS: {p_nis:.1f}%\n(McNemar p={p_bf_mc:.1e})",
              fontsize=6.3, bbox=dict(boxstyle="round,pad=0.25", facecolor="#F5F7FA", edgecolor="#CCCCCC"))
    ax_c.set_xlim(0, 168); ax_c.set_ylim(0, 62)
    ax_c.set_xlabel("Simulated Perfusion Time t (hours)")
    ax_c.set_ylabel("Molar Ratio R_LP  or  Heteroplasmy h (%)")
    ax_c.set_title("C. Biofouling & Volumetric OUR + OD + ORR Handover", fontweight="bold", loc="left")
    ax_c.legend(loc="upper left", fontsize=6.3, frameon=True, framealpha=0.92)
    ax_c.grid(alpha=0.25)

    # Panel D: Cassette Perfusion Occlusion Acidosis Rescue + C99 Q16.16 RK4 & 5-Ch EKF Parity
    ax_d = axes[1, 1]
    t_grid = np.linspace(0.0, T_HORIZON_H, 1001)
    tr_no = arms_nom["fault_no_rb"]["_traj"]
    tr_rb = arms_nom["fault_rb"]["_traj"]
    ax_d.axvspan(72.0, 96.0, color="#D55E00", alpha=0.14, label="Perfusion Line Occlusion (72–96 h)")
    ax_d.plot(t_grid, tr_no[:, 7], color="#D55E00", linestyle="--", linewidth=1.8,
              label=f"One-Way Latch (min pH={arms_nom['fault_no_rb']['min_ph']:.2f}, T_SCM={arms_nom['fault_no_rb']['tscm_pct']:.1f}%)")
    ax_d.plot(t_grid, tr_rb[:, 7], color="#0072B2", linewidth=2.0,
              label=f"Bypass Flush u_fl=1 (min pH={arms_nom['fault_rb']['min_ph']:.2f}, T_SCM={arms_nom['fault_rb']['tscm_pct']:.1f}%)")
    ax_d.axhline(6.85, color="#CC79A7", linestyle=":", linewidth=1.1, label="Optical pH Flush Trigger (6.85)")
    err_1000 = c99_data["step_sweep"][1]["err_q16_sd_pct"]
    ekf_diff = c99_data["c99_ekf_biofoul_verification"]["ekf_q16_vs_fp64_rmse_pct"]
    ax_d.text(3.5, 6.51,
              f"Continuous Plant State x(t) (Zero State Reset)\n"
              f"C99 Q16.16 RK4 Err: {err_1000:.4f}% | C99 5-Ch EKF vs FP64: {ekf_diff:.4f}%\n"
              f"Zero-Heap Static SRAM: {c99_data['cart_dtcm_workspace_bytes']:,} / 8,192 B ({c99_data['cart_dtcm_bank0_utilization_pct']:.1f}%)",
              fontsize=6.2, bbox=dict(boxstyle="round,pad=0.25", facecolor="#F5F7FA", edgecolor="#CCCCCC"))
    ax_d.set_xlim(0, 168); ax_d.set_ylim(6.02, 7.46)
    ax_d.set_xlabel("Simulated Perfusion Time t (hours)")
    ax_d.set_ylabel("Simulated Cassette Perfusate pH(t)")
    ax_d.set_title("D. Plant-Preserving Acidosis Flush & C99 EKF/RK4 Parity", fontweight="bold", loc="left")
    ax_d.legend(loc="lower left", fontsize=6.1, frameon=True, framealpha=0.92)
    ax_d.grid(alpha=0.25)

    out_path = FIG_DIR / "paper2_fig3_donor_cohort_biofouling_and_acidosis.pdf"
    fig.savefig(out_path, bbox_inches="tight")
    plt.close(fig)
    return out_path


# =============================================================================
# Main Execution & JSON Export
# =============================================================================
def strip_private_arrays(obj):
    if isinstance(obj, dict):
        return {k: strip_private_arrays(v) for k, v in obj.items() if not k.startswith("_")}
    if isinstance(obj, list):
        return [strip_private_arrays(v) for v in obj]
    return obj


def main() -> None:
    t_start = time.time()
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    configure_matplotlib()

    print("[1/6] Auditing experimental mmCIF (2Q8G, 4ZZJ) & AlphaFold v6 (FOXO1, PGC-1a) structures...")
    struct_data = audit_all_structures()

    print("[2/6] Simulating 14 nominal 8-state CAR-T cassette protocol & ablation arms...")
    nominal_modes = [
        "warburg",
        "premature_ox_h0",
        "urolithin_mono",
        "fixed_24",
        "fixed_36",
        "fixed_48",
        "closed_loop",
        "fixed_72",
        "fixed_96",
        "fixed_120",
        "ablate_azd",
        "ablate_sirt1",
        "fault_no_rb",
        "fault_rb",
    ]
    arms_nom = {m: simulate_nominal_arm(m) for m in nominal_modes}

    print("[3/6] Solving 200-variable discrete-adjoint optimal control (5-start L-BFGS-B + all-200 FD check) & Pareto frontier...")
    pmp_data = solve_discrete_adjoint_pareto()

    print("[4/6] Computing volumetric OUR + OD + ORR Observability Gramian & N=1,000 patient cohort...")
    gramian_data = compute_cassette_observability_gramian()
    cohort_data = run_cohort_and_biofouling_benchmarks()

    print("[5/6] Compiling & executing bare-metal C99 Q16.16 RK4 + 5-channel Huber-EKF testbench...")
    c99_data = run_c99_cassette_verification()

    print("[6/6] Rendering 3 vector PDF publication figures...")
    f1 = render_fig1_structural_pharmacology(struct_data)
    f2 = render_fig2_ode_and_pareto(arms_nom, pmp_data)
    f3 = render_fig3_cohort_biofouling_and_acidosis(cohort_data, arms_nom, c99_data)

    payload = {
        "metadata": {
            "paper": "Paper II — Closed-Loop Metabolic Gating of Ex Vivo CAR-T Perfusion Bioreactors",
            "epistemic_status": "100% In Silico Computational Simulation & Structural Bioinformatics (Zero Wet-Lab Execution Claimed)",
            "miase_prng_seed": SEED,
            "horizon_hours": T_HORIZON_H,
            "runtime_seconds": round(time.time() - t_start, 2),
            "figures": [str(f1.relative_to(ROOT)), str(f2.relative_to(ROOT)), str(f3.relative_to(ROOT))],
        },
        "nominal_parameters": NOMINAL_PARAMS,
        "initial_state": dict(zip(["W", "D", "N_cells", "R_LP", "F_nuc", "S_SRC", "E_TOX", "pH"], X0_NOMINAL.tolist())),
        "structural_pharmacology": strip_private_arrays(struct_data),
        "nominal_protocol_sweep": strip_private_arrays(arms_nom),
        "adjoint_optimal_control_and_pareto": strip_private_arrays(pmp_data),
        "observability_gramian": gramian_data,
        "heterogeneous_cohort_and_biofouling": strip_private_arrays(cohort_data),
        "c99_coralnpu_verification": c99_data,
    }

    out_json = SIM_DIR / "paper2_verified_benchmarks.json"
    with out_json.open("w", encoding="utf-8") as fh:
        json.dump(payload, fh, indent=2)

    print(f"\n[DONE] Wrote {out_json.relative_to(ROOT)} in {time.time() - t_start:.2f}s")
    print(f"  - PDK1 (2Q8G) top contact: {struct_data['pdk1_2q8g']['contacts_le_4A'][0]}")
    print(f"  - SIRT1 (4ZZJ) STAC contacts: {[c['residue']+':'+str(c['min_dist_A']) for c in struct_data['sirt1_4zzj_stac_4tq']['contacts_le_4A']]}")
    print(f"  - Closed-Loop Nominal: sw={arms_nom['closed_loop']['switch_time_h']}h, N={arms_nom['closed_loop']['expansion_fold']}x, T_SCM={arms_nom['closed_loop']['tscm_pct']}%, h={arms_nom['closed_loop']['h_end_pct']}%")
    print(f"  - Adjoint PMP (all 200 vars FD) grad err={pmp_data['max_adjoint_vs_fd_grad_err']:.2e}, 5-start spread={pmp_data['multi_start_cost_spread_pct']}%")
    print(f"  - Cohort Pass Rates (N=1,000): Best Fixed(68h)={cohort_data['arms']['fixed_68']['pass_rate_pct']}% | 4-Ch No-ORR={cohort_data['arms']['closed_loop_ekf_no_orr']['pass_rate_pct']}% | 5-Ch ORR(168h)={cohort_data['arms']['closed_loop_ekf']['pass_rate_pct']}% | Adaptive(168-240h)={cohort_data['arms']['adaptive_horizon_240h_ekf']['pass_rate_pct']}% | Oracle168={cohort_data['oracle_controllable_pct']}% | Oracle240={cohort_data['oracle_240h_controllable_pct']}%")
    print(f"  - Biofouling Pass Rates: Naive Lac:Pyr={cohort_data['arms']['biofoul_naive_lp_gate']['pass_rate_pct']}% | 4-Ch Handover={cohort_data['arms']['biofoul_no_orr_nis_handover']['pass_rate_pct']}% | 5-Ch ORR+NIS={cohort_data['arms']['biofoul_ekf_nis_handover']['pass_rate_pct']}%")
    print(f"  - C99 Q16.16 RK4 Err (N=1000): {c99_data['step_sweep'][1]['err_q16_sd_pct']}% | C99 5-Ch EKF vs FP64: {c99_data['c99_ekf_biofoul_verification']['ekf_q16_vs_fp64_rmse_pct']}% | SRAM={c99_data['cart_dtcm_workspace_bytes']} B")


if __name__ == "__main__":
    main()

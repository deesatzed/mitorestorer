---
name: bioenergetic-precision-protocol
description: >-
  End-to-end precision bioenergetic discovery, biophysical ODE simulation, 3D
  atomic structural validation, and biomarker-gated clinical protocol generator
  grounded in Picard & Murugan's Energy Resistance Principle (EP = eR * f^2)
  and the Microbiota-Gut-Mitochondria Axis. Use when the user provides a target
  disease/syndrome (e.g., metabolic syndrome, neurodegenerative disease,
  refractory epilepsy, treatment-resistant depression, Long COVID/ME-CFS)
  and/or host genetic variants (SNPs) and wants to simulate mitochondrial
  heteroplasmy clearance vs. biogenesis, validate cofactor/drug rescue in 3D,
  or design a chrono-phased synbiotic and nutraceutical protocol.
---

# Bioenergetic Precision Protocol & ODE Simulator (`bioenergetic-precision-protocol`)

## Overview
Orchestrates a 4-stage precision medicine workflow that couples **Picard & Murugan's Energy Resistance Principle ($\text{EP} = \text{éR} \cdot f^2$)** with the **Microbiota–Gut–Mitochondria Axis**, **Host Genomic Variant Parameterization**, **Sub-Hourly Biophysical ODE Simulation (`dt = 0.02`d)**, and **3D Atomic Structural Validation**:

1. **Stage A — Empirical Multi-Omic Parameterization:** Resolves host genetic variants (`GRCh38` coordinates, canonical MANE Select transcripts, `AlphaMissense`, `Mutfunc ddG`, `LOEUF`, `Haplo/Triplo` dosage sensitivity, `gnomAD v4` allele frequencies, and `GTEx`/`OpenTargets` tissue profiles) into quantitative ODE rate multipliers.
2. **Stage B — Coupled Biophysical ODE & Sensitivity Simulation:** Runs the sub-hourly (`dt = 0.02` days) coupled differential equation engine (`scripts/bcr_pipeline.py simulate`) modeling wild-type ($W$) vs. deletion-mutant ($D$) mtDNA heteroplasmy, hepatic glycogen / `PINK1` auto-phosphorylation lag ($\tau_{\text{lag}}$), gut $\text{H}_2\text{S}$ Complex IV inhibition vs. butyrate/Urolithin A/bacterial *PncA* $\text{NAD}^+$ salvage, `TOMM40` import impedance, synaptic GABA/Glutamate excitability ($\text{EI}_{\text{ratio}}$), microglial `IDO1`/Quinolinic Acid (`QUIN`), and systemic `GDF15` gating.
3. **Stage C — 3D Atomic Structural Validation:** Downloads experimental X-ray / Cryo-EM structures via `pdb-database` (or falls back to `alphafold-database-fetch-and-analyze`), computes exact inter-atomic contact distances (`< 4.5 Å`) via `scripts/bcr_pipeline.py inspect-cif`, and loads + annotates the structures in the interactive [Biomolecular Structure Viewer](sidecar://science/biomolecular_structure_viewer/ui.html).
4. **Stage D — Deliverable Synthesis:** Generates an interactive HTML simulator (`generative_ui`), a biomarker-gated clinical trial protocol, and a clinician/patient field manual & lab gating scorecard (`scripts/bcr_pipeline.py verify-gate`).

---

## Dependencies
This skill coordinates and references the following installed skills (never reimplement their APIs directly):
- **`uv`**: Required to execute `uv run` for all helper scripts.
- **`ensembl-database`**: Resolves gene/variant IDs, MANE Select transcripts, and VEP pathogenicity/structural scores (`AlphaMissense`, `Mutfunc ddG`, `Dosage Sensitivity`).
- **`gnomad-database`**: Retrieves population allele frequencies (`AF`) and gene constraint metrics (`LOEUF`, `pLI`).
- **`opentargets-database`**: Queries target-disease associations and mitokine (`GDF15`, `FGF21`) clinical phenotypes.
- **`gtex-database`**: Retrieves organ-specific baseline expression (`TPM`) to scale tissue turnover ($\tau_{\text{tissue}}$).
- **`pdb-database`**: Searches and downloads experimental `.cif` 3D structures (`fetch_pdb_metadata.py`, `download_coordinate_files.py`).
- **`alphafold-database-fetch-and-analyze`**: Automatic fallback when an experimental structure is unavailable in `pdb-database`.
- **`biomolecular-structure-viewer`**: Renders `.cif` files and attaches 3D residue/ligand pins (`load_structure`, `annotate_region`, `set_representation`).
- **`generative_ui`**: Provides design guidelines for building the standalone interactive HTML ODE simulator.

---

## Quick Start

```bash
# 1. Run a 90-day coupled ODE simulation across all 4 clinical cohorts and write telemetry to JSON
uv run --no-project /Users/o2satz/.gemini/config/skills/bioenergetic-precision-protocol/scripts/bcr_pipeline.py simulate \
  --cohort all \
  --pink1-class Class_I_ATP_Pocket \
  --days 90 \
  --output bcr_sim_telemetry.json

# 2. Inspect a decompressed .cif file for chains, ligands, and atomic contacts (< 4.5 A) around a residue/ligand
uv run --no-project /Users/o2satz/.gemini/config/skills/bioenergetic-precision-protocol/scripts/bcr_pipeline.py inspect-cif \
  --cif 4ZZJ/pdb_00004zzj.cif \
  --target-chain A \
  --target-res 703 \
  --max-dist 4.5 \
  --output sirt1_stac_contacts.json

# 3. Evaluate Day-38 or Day-54 patient lab telemetry for Stage 1 -> Stage 2 transition gating
uv run --no-project /Users/o2satz/.gemini/config/skills/bioenergetic-precision-protocol/scripts/bcr_pipeline.py verify-gate \
  --gdf15-day0 1420 \
  --gdf15-current 840 \
  --day 38 \
  --lean-mass-delta-pct -0.9 \
  --bhb-mmol 1.8 \
  --output gate_decision_day38.json
```

---

## Utility Scripts

All programmatic work is handled by `/Users/o2satz/.gemini/config/skills/bioenergetic-precision-protocol/scripts/bcr_pipeline.py`. All subcommands write JSON output to `--output` and print a brief status message to stdout.

### 1. `simulate`
Integrates the sub-hourly (`dt = 0.02`d) coupled ODE system across peripheral metabolic, neurodegenerative (`APOE4/TOMM40`), refractory epilepsy (`PV+` interneuron $\text{EI}_{\text{ratio}}$ & Day-30 $\text{H}_2\text{S}$ escape), and treatment-resistant depression (`IDO1/QUIN` & `GDF15–GFRAL` conservation shutdown) cohorts.
- Required flags: `--cohort` (`metabolic_muscle`, `neuro_apoe4`, `refractory_epilepsy`, `treatment_resistant_mdd`, or `all`), `--days` (integer, e.g., `90`), `--output` (filepath).
- Optional flags: `--pink1-class` (`WT`, `Class_I_ATP_Pocket`, `Class_II_Ins3_Ub_Cradle`).

### 2. `inspect-cif`
Parses a local decompressed `.cif` file and extracts author chains (`auth_asym_id`), modified residues (`SEP`, `TPO`, `ALY`), non-polymer cofactors/ligands, and all inter-atomic contacts within `--max-dist` Angstroms of a specified residue or ligand.
- Required flags: `--cif` (path to `.cif`), `--max-dist` (float, e.g. `4.5`), `--output` (filepath).
- Optional flags: `--target-chain` (author chain ID), `--target-res` (author residue number or `auth_comp_id`).

### 3. `verify-gate`
Applies the Two-Stage Biomarker-Gated "Traffic Light" rules to patient lab telemetry (`GDF15`, capillary $\beta$-OHB, lean mass change) to determine whether to transition from Stage 1 (Lag-Compressed 18:6 Purge) to Stage 2 (Weekly `2.25d Purge : 4.75d Build` Expansion).
- Required flags: `--gdf15-day0`, `--gdf15-current`, `--day`, `--lean-mass-delta-pct`, `--output`.
- Optional flags: `--bhb-mmol` (default `1.5`).

---

## Workflow

1. **Read Reference Equations & Calibration Table:** Read `/Users/o2satz/.gemini/config/skills/bioenergetic-precision-protocol/references/erp_ode_reference.md` to review the 5 core kinetic discoveries and empirical SNP multipliers.
2. **Stage A (Multi-Omic Variant & Disease Parameterization):**
   - If the user supplies a disease/syndrome and/or SNPs, use `opentargets-database`, `ensembl-database`, `gnomad-database`, and `gtex-database` to fetch empirical variant consequences and tissue expression, or map them to the pre-calibrated archetypes in `references/erp_ode_reference.md`.
3. **Stage B (Sub-Hourly Coupled ODE Simulation):**
   - Run `uv run .../scripts/bcr_pipeline.py simulate` with `--output` saved in the conversation's scratch/data directory. Inspect the JSON output for Day-90 heteroplasmy ($h\%$), ATP capacity, GDF15 trajectory, switch day, and any `"safety_warnings"`.
4. **Stage C (3D Atomic Structural Proof):**
   - Read `pdb-database` and `biomolecular-structure-viewer` skills. Ensure `.licenses/pdb_database_LICENSE.txt` notification is satisfied.
   - Download target `.cif.gz` coordinate files into the conversation's artifact directory via `pdb-database` (or `alphafold-database-fetch-and-analyze` if no experimental PDB entry exists), decompress with `gunzip -f`, and run `uv run .../scripts/bcr_pipeline.py inspect-cif` to verify author residue numbers (`auth_asym_id`, `auth_seq_id`) and measure exact inter-atomic distances (`< 4.5 Å`).
   - Load and annotate the structures in the viewer using `load_structure`, `annotate_region`, and `set_representation`, and **always** include `[Biomolecular Structure Viewer](sidecar://science/biomolecular_structure_viewer/ui.html)` in the chat response.
5. **Stage D (Protocol & Interactive Simulator Synthesis):**
   - Generate or update the interactive HTML simulator, clinical trial protocol, and clinician/patient field manual artifacts.

---

## Rate Limiting
- `scripts/bcr_pipeline.py` performs purely local mathematical ODE integration and local `.cif` parsing (zero network requests).
- All external genomic and structural database queries must go through the referenced skills (`ensembl-database`, `gnomad-database`, `opentargets-database`, `gtex-database`, `pdb-database`, `alphafold-database-fetch-and-analyze`), which enforce cross-process file-lock rate limiting via `polite-http`.

---

## Common Mistakes
1. **Passing `label_seq_id` instead of `auth_seq_id` to the 3D Viewer:** `biomolecular-structure-viewer` strictly requires author numbering (`auth_asym_id`, `auth_seq_id`). Always run `bcr_pipeline.py inspect-cif` first to verify exact author IDs.
2. **Co-Administering $\text{NAD}^+$ Precursors During Catabolic Purge Windows (Discovery 1):** Never schedule NR/NMN or high-leucine protein during Stage 1 fasting/Purge windows while heteroplasmy $h > 8\%$, as shorter deletion-mutant mtDNA ($D$) replicates ~32% faster than wild-type $W$.
3. **Treating Class II `PINK1` Ubiquitin-Interface Mutants with Kinetin (`KTP`) Monotherapy (`PDB: 6EQI`):** Variants disrupting the `Insertion-3` Ubiquitin cradle (`residues 265–285`) or `SEP202–Arg282` clamp cannot dock `Ubiquitin Ser65` even when `KTP` occupies the ATP pocket; they require `Trehalose + Spermidine` (`TFEB` macro-autophagy bypass) and nocturnal N3 astrocytic transmitophagy.

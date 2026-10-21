# Phase 2 Skill Design Plan: `bioenergetic-precision-protocol`

## 1. Skill Name & YAML Frontmatter Description

```yaml
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
```
* **Name Length:** 31 characters ($\le 64$ limit, lowercase + hyphens).
* **Description Length:** 578 characters ($\le 1024$ limit).

---

## 2. Directory Structure & Global Install Location

Installed globally for Antigravity 2.0 at `/Users/o2satz/.gemini/config/skills/bioenergetic-precision-protocol/` so it is automatically discovered in all future conversations:

```text
/Users/o2satz/.gemini/config/skills/bioenergetic-precision-protocol/
├── SKILL.md                        # Main workflow orchestration & rules
├── scripts/
│   └── bcr_pipeline.py             # Pure-stdlib CLI engine (simulate, inspect-cif, verify-gate)
└── references/
    └── erp_ode_reference.md        # Equations, 5 kinetic discoveries, & SNP-to-ODE calibration table
```

---

## 3. Existing Skills Referenced (Rule 1 Compliance)

The new skill delegates all external database queries and 3D rendering to existing installed skills rather than reimplementing them:

| Referenced Skill | Role in Workflow |
| :--- | :--- |
| **`uv`** | Ensures `uv run` is available to execute `scripts/bcr_pipeline.py` and upstream database scripts. |
| **`ensembl-database`** | Resolves rsIDs and genomic coordinates (`GRCh38`), canonical MANE Select transcripts, `AlphaMissense`, `SIFT`, `PolyPhen`, `Mutfunc ddG`, and `Dosage Sensitivity (Haplo/Triplo)`. |
| **`gnomad-database`** | Retrieves global/population allele frequencies (`AF`), homozygote counts, and gene loss-of-function constraint (`LOEUF`, `pLI`). |
| **`opentargets-database`** | Retrieves target-disease association scores, tractability, and systemic mitokine (`GDF15`, `FGF21`) clinical phenotypes. |
| **`gtex-database`** | Retrieves tissue-specific gene expression (`TPM`) across brain, muscle, liver, and gut to calibrate organ turnover ($\tau_{\text{tissue}}$). |
| **`pdb-database`** | Primary source for experimental X-ray / Cryo-EM 3D structures (`fetch_pdb_metadata.py` and `download_coordinate_files.py`). |
| **`alphafold-database-fetch-and-analyze`** | Automatic fallback when a target human protein lacks an experimental structure in `pdb-database`. |
| **`biomolecular-structure-viewer`** | Renders downloaded `.cif` files and attaches 3D residue/ligand annotations (`load_structure`, `annotate_region`, `set_representation`) in the [Biomolecular Structure Viewer](sidecar://science/biomolecular_structure_viewer/ui.html). |
| **`generative_ui`** | Guides the styling and layout of the standalone interactive HTML ODE simulator artifact (`bcr_simulator.html`). |

---

## 4. CLI Script Specification (`scripts/bcr_pipeline.py`)

Follows **Rule 3 (CLI Script Pattern)** and **Rule 4 (Default to File Output)** using only Python standard library modules (`argparse`, `json`, `math`, `os`, `sys`). Every subcommand requires `--output <filepath>` and writes structured JSON to disk while printing only a brief status confirmation to stdout.

### Subcommand 1: `simulate`
Runs the sub-hourly (`dt = 0.02` days = `0.48` hours) coupled Bioenergetic / Microbiome / Transmitophagy ODE solver across regimens and cohorts.
* **Arguments:**
  * `--cohort` (required): `metabolic_muscle` | `neuro_apoe4` | `refractory_epilepsy` | `treatment_resistant_mdd` | `all`
  * `--pink1-class` (optional, default `Class_I_ATP_Pocket`): `WT` | `Class_I_ATP_Pocket` | `Class_II_Ins3_Ub_Cradle`
  * `--days` (required): Simulation horizon in days (e.g., `90`)
  * `--output` (required): Output JSON file path

### Subcommand 2: `inspect-cif`
Parses a local decompressed `.cif` / `.mmcif` coordinate file, lists polymer chains (`auth_asym_id`), modified residues (`SEP`, `TPO`, `ALY`, etc.), non-polymer ligands (`FAD`, `CNA`, `4TQ`, `PC1`, `ZN`), and computes exact inter-atomic contact distances (`< --max-dist` Å) around a target residue or ligand for `annotate_region`.
* **Arguments:**
  * `--cif` (required): Path to local `.cif` file
  * `--target-chain` (optional): Author chain ID (e.g., `A`, `C`)
  * `--target-res` (optional): Author residue number or ligand `auth_comp_id` (e.g., `222`, `202`, `CNA`, `4TQ`)
  * `--max-dist` (required): Distance cutoff in Angstroms (e.g., `4.5`)
  * `--output` (required): Output JSON file path

### Subcommand 3: `verify-gate`
Evaluates Day-38 or Day-54 patient laboratory telemetry (`GDF15`, capillary `BHB`, `Lactate:Pyruvate`, lean mass delta) and outputs the Green/Yellow/Red Stage 1 $\rightarrow$ Stage 2 gating decision.
* **Arguments:**
  * `--gdf15-day0` (required): Baseline fasting plasma GDF15 ($\text{pg/mL}$)
  * `--gdf15-current` (required): Current fasting plasma GDF15 ($\text{pg/mL}$)
  * `--day` (required): Current protocol day (e.g., `38` or `54`)
  * `--lean-mass-delta-pct` (required): Percent change in lean mass from Day 0 (e.g., `-0.8`)
  * `--bhb-mmol` (optional, default `1.5`): Morning Hour-15 capillary $\beta$-OHB ($\text{mmol/L}$)
  * `--output` (required): Output JSON file path

---

## 5. Rate Limiting Strategy (Rule 2 Compliance)

* All external biological API calls (`Ensembl`, `gnomAD`, `OpenTargets`, `GTEx`, `RCSB PDB`, `AlphaFold DB`) are executed exclusively through the existing referenced skills (`ensembl-database`, `gnomad-database`, `opentargets-database`, `gtex-database`, `pdb-database`, `alphafold-database-fetch-and-analyze`), which already enforce cross-process file-lock rate limiting via `polite-http`.
* `scripts/bcr_pipeline.py` performs purely local mathematical ODE integration and local `.cif` file parsing (zero direct external HTTP calls), eliminating any unguarded network rate-limit risk.

---

## 6. Per-Step Error Handling Strategy

1. **Stage A (Genomic & Multi-Omic Lookup):** If a user provides only a disease name without specific SNPs, query `opentargets-database` first to identify top bioenergetic/mitochondrial genetic targets for that disease, then resolve coordinates via `ensembl-database`. If an API is unreachable, surface the exact error and allow fallback to the pre-calibrated reference table in `references/erp_ode_reference.md`.
2. **Stage B (ODE Simulation):** If a simulated regimen triggers cellular mass atrophy ($W + D < 0.75$) or Day-30 sulfidogenic escape ($\text{H}_2\text{S} > 1.2$), `bcr_pipeline.py simulate` flags a structured `"safety_warnings"` array in the output JSON so the agent automatically adjusts the Purge:Build ratio or adds the Neuro-BCR / Sulfide-Guard stack.
3. **Stage C (3D Structural Validation):**
   * Always check for `.licenses/pdb_database_LICENSE.txt` and notify the user of https://www.rcsb.org/pages/usage-policy per `pdb-database`.
   * If `pdb-database` returns no experimental structure for a target protein, automatically fall back to `alphafold-database-fetch-and-analyze` (resolving the UniProt accession via `uniprot-database` / `ensembl-database`) and inform the user that an AlphaFold predicted model is being displayed.
   * Always decompress `.cif.gz` files before calling `load_structure`, verify author numbering (`auth_asym_id`, `auth_seq_id`) via `bcr_pipeline.py inspect-cif`, and include `[Biomolecular Structure Viewer](sidecar://science/biomolecular_structure_viewer/ui.html)` in the chat response.

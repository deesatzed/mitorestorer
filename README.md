# MitoRestorer: Precision Bioenergetic Circuit Reset (BCR) Platform

**MitoRestorer** is an end-to-end biophysical discovery, sub-hourly ordinary differential equation (ODE) simulation, 3D atomic structural validation, and biomarker-gated clinical protocol platform grounded in **Picard & Murugan’s Energy Resistance Principle ($\text{EP} = \text{éR} \cdot f^2$)**, the **Microbiota–Gut–Mitochondria Axis**, and **Host Genetic Variant Parameterization** (`Ensembl VEP`, `AlphaMissense`, `Mutfunc`, `gnomAD v4`, `OpenTargets`, and `GTEx v8`).

### Key Entry Points
* **Interactive Landing Page:** [`index.html`](index.html)
* **Interactive 90-Day Web Simulator (with 4-Step Guided Wizard & Lab Gate Calculator):** [`simulator/bcr_simulator.html`](simulator/bcr_simulator.html)
* **Beginner / College-Freshman Step-by-Step User Guide:** [`docs/user_walkthrough.md`](docs/user_walkthrough.md)
* **Clinician & Patient Daily Field Manual:** [`docs/bcr_field_manual.md`](docs/bcr_field_manual.md)
* **Full 5-Cohort Clinical Protocol Specification:** [`docs/bcr_clinical_protocol.md`](docs/bcr_clinical_protocol.md)
* **End-to-End Executable Sample Patient Script:** [`examples/run_sample_patient_walkthrough.py`](examples/run_sample_patient_walkthrough.py)

---

## Step-by-Step User Walkthrough (Beginner / College-Freshman Guide)

Think of **mitochondria** as thousands of tiny rechargeable batteries inside every cell. When a percentage of those batteries suffer DNA damage ($h\%$ / **heteroplasmy**), they resist electron flow (**Energy Resistance**, $\text{éR}$). Under Picard & Murugan’s **Energy Resistance Principle ($\text{EP} = \text{éR} \cdot f^2$)**, forcing high metabolic flux ($f$) through a high-resistance circuit causes cellular stress ($\text{EP}$) to spike with the **square** of the flux, triggering a systemic blood stress hormone called **GDF15** (the cell’s "smoke alarm").

Here is how to run a complete patient case from start to finish in **4 steps** (either inside [`simulator/bcr_simulator.html`](simulator/bcr_simulator.html) by clicking **`▶ Load & Run Sample Patient (Step 1)`**, or in the terminal via `python3 examples/run_sample_patient_walkthrough.py`):

### Step 1 of 4: Pick a Clinical Scenario & Inspect Genetic Bottlenecks
* Select one of the **5 pre-calibrated clinical cohorts** (e.g., **Cohort E: Early Parkinson’s Disease**, sample patient `SAMPLE-PD-001` starting at `58.0%` damaged batteries and `GDF15 = 1,950 pg/mL`) or customize the **Empirical SNP Panel**:
  * **`PPARGC1A rs8192678`** (`[Battery Factory Builder]`): Controls how fast the cell builds new mitochondria (`0.58x–1.00x`).
  * **`FUT2 rs601338`** (`[Gut Microbiome Feeder]`): `A/A Non-Secretors` (**22.5% of people**) cannot convert dietary pomegranate/berries into the battery-cleaning postbiotic **Urolithin A**.
  * **`MTHFR rs1801133`** (`[Methylation & Creatine Supply]`): Drains `SAMe` and Creatine unless rescued by **Riboflavin B2 (`PDB: 6FCX`) + Creatine**.
  * **`APOE4 / TOMM40 rs429358`** (`[Mitochondrial Import Door]`): Clogs the outer-membrane `TOMM40` pore (`PDB: 7CK6`) unless unblocked by **TUDCA + Citicoline/LPC-DHA**.
  * **`PINK1 / PRKN`** (`[Damaged Battery Trash Collector]`): Tags damaged mitochondria for removal (`PDB: 6EQI`).
  * **`GBA1 rs76763715` & `LRRK2 rs34637584`** (`[Lysosomal Recycling Gate]`): Traps glucocerebrosidase in the ER (`PDB: 2NSX`) and blocks autophagosome–lysosome fusion in Parkinson’s disease.

### Step 2 of 4: See Why Standard Treatments & "All-at-Once" Stacks Fail
* Click **`Next Step →`** (or select **Regimen 4: Unphased "All-in-One" Stack**):
  * **Discovery 1 (The Simultaneous Biogenesis Trap):** Taking "cleaning" supplements (`Urolithin A`) and "building" supplements (`NAD+ / NR` + high protein) **at the same time** keeps the cell's growth switch (`mTORC1`) ON. That locks the trash collector (`ULK1`) shut while telling `POLG` DNA polymerase to copy existing batteries. Because broken mtDNA loops ($D$) are shorter, they copy **32% faster** than healthy ones—cloning broken batteries to **>84%**!
  * **Discovery 6 (The `GBA1/LRRK2` Autophagosome Traffic Jam):** In Parkinson’s disease, stimulating mitophagy (`Urolithin A + Kinetin`) *before* unblocking `GBA1` lysosomes causes wrapped-up trash bags to pile up (`Auto_Jam = 1.29x`), worsening **Dopamine-Quinone (`DA-Q`)** and $\alpha$-synuclein clumping (`1.63x`).

### Step 3 of 4: Activate the 3-Stage Phased Reset & 3D Cofactor Rescue
* Click **`Next Step →`** (switching to **Regimen 5: Two-Stage Biomarker-Gated BCR**):
  1. **Stage 0 (Days 1–14 — Unblock the Doors & Incinerators First):** High-Dose **Ambroxol** (`PDB: 2NSX`) + **Spermidine** (`ATP13A2` substrate) + **TUDCA + Citicoline** (`PDB: 7CK6`) + **Molybdenum/Polyphenols/Tributyrin** (blocks gut $\text{H}_2\text{S}$ toxin).
  2. **Stage 1 (Days 15–38 — Lag-Compressed Daily 18:6 Cleaning Purge):** Standard 16:8 fasting wastes `13.5 hours` burning liver glycogen before mitophagy starts. Combining a **morning Zone-2 walk + C8-MCT Oil + Exogenous D-BHB** with **Direct Urolithin A (`1,000 mg`) + Kinetin (`PDB: 6EQI`)** compresses that lag to **`5–6 hours`**, clearing damaged mitochondria from **`58.0% → 3.7%`** while holding `NR/NMN` OFF.
  3. **Stage 2 (Post-Gate → Day 90 — Weekly `2.25d Purge : 4.75d Build`):** Once broken batteries are <8%, activates **NR (`500–1,000 mg`) + *PncA* Synbiotic + Resveratrol (`SIRT1 PDB: 4ZZJ`) + Riboflavin B2 (`PDB: 6FCX`) + Creatine (`5 g`)** for 4.75 days/week to clone healthy wild-type ($W$) mitochondria (`ATP = 1.20x–1.46x`).

### Step 4 of 4: Verify the Day-38 Blood Test Gate & Daily Action Plan
* Click **`Next Step →`** (opening the **`Action Plan & Lab Gate`** tab):
  * Enter the patient's **Day-0 (`1,950 pg/mL`)** and **Day-38 (`510 pg/mL`, `-73.8%` drop)** fasting plasma **GDF15**, **Morning Ketones (`1.8 mmol/L`)**, and **Lean Mass Change (`-0.6%`)**.
  * The calculator outputs **`🟢 GREEN LIGHT: Flip to Stage 2 Build`** (since GDF15 dropped $\ge 35\%$ with $<2.5\%$ muscle loss) and renders the patient's personalized **Stage 0 / Stage 1 / Stage 2 Daily Checklist**.

---

## 5 Clinical Cohorts & 90-Day Simulation Summary

| Cohort & Indication | Genetic & Circuit Bottleneck | Standard / Unphased Care (Day 90) | Phased BCR Protocol (Day 90) |
| :--- | :--- | :--- | :--- |
| **Cohort A: Metabolic / Sarcopenia** | `PINK1` Class I + `FUT2 rs601338 A/A` Non-Secretor | 16:8 Fast: $h = 69.4\%$, $\text{ATP} = 0.50\times$<br>Unphased: $h = 91.4\%$ | **WT:** $h = 3.3\%$, **$\text{ATP} = 1.20\times$**<br>**Risk:** $h = 20.0\%$, **$\text{ATP} = 1.02\times$** |
| **Cohort B: Prodromal *APOE4* AD** | `APOE4 rs429358` (`1B68`) + `TOMM40` stall (`7CK6`) | Peripheral BCR only: $h = 94.2\%$, $\text{ATP} = 0.17\times$ | **`Neuro-BCR`:** $h = 3.3\%$, **$\text{ATP} = 1.19\times$**, Mass $= 1.06\times$ |
| **Cohort C: Refractory Epilepsy** | Fast-spiking $\text{PV}^+$ interneurons ($f^2 = 3.06\times$) + Keto $\text{H}_2\text{S}$ | Classic KD: Day-33 $\text{H}_2\text{S} = 1.97$, $\text{E/I} = 2.26$ | **`Epi-BCR`:** $h = 0.5\%$, $\text{H}_2\text{S} = 0.11$, **$\text{E/I} = 0.65$ ($96\%$ seizure red.)** |
| **Cohort D: Treatment-Resistant MDD** | Microglial `IDO1` $\rightarrow$ `QUIN` + Brainstem `GDF15–GFRAL` lock | Ketamine Alone: Spines $0.57\times$ (d7) $\rightarrow 0.45\times$ (d90) | **`Affective-BCR`:** $h = 0.6\%$, $\text{QUIN} \downarrow 87\%$, **Spines $= 1.24\times$** |
| **Cohort E: Parkinson's Disease (`SNc`)** | `GBA1 rs76763715` (`2NSX`) + `LRRK2 rs34637584` + `DA-Q` | Mitophagy w/o Ambroxol: $\text{Jam} = 1.29\times$, $\text{SNCA} = 1.63\times$ | **`PD-BCR`:** $h = 3.7\%$, $\text{Jam} = 0.08\times$, **$\text{TH Tone} = 1.25\times$ ($95\%$ rescue)** |

---

## Repository Structure

```text
mitorestorer/
├── index.html                       # Interactive Project Landing Page & 4-Step Walkthrough
├── README.md                        # Overview, Step-by-Step Guide, & Quick Start
├── docs/
│   ├── user_walkthrough.md          # Beginner / College-Freshman Step-by-Step User Walkthrough
│   ├── bcr_clinical_protocol.md     # Complete 5-Cohort Clinical Protocol & Phase IIa Trial Design
│   ├── bcr_field_manual.md          # Clinician & Patient Daily Checklist + Day 38/54 GDF15 Lab Worksheet
│   ├── research_roadmap.md          # Multi-Omic & Structural Validation Roadmap
│   └── skill_design_plan.md         # Architecture of the Reusable Agent Skill
├── simulator/
│   ├── bcr_simulator.html           # Interactive 90-Day HTML5/Canvas ODE Simulator (Cohorts A–E)
│   ├── sim_bioenergetic_reset.py    # Core Python Biophysical ODE Solver across all 5 Cohorts
│   ├── analyze_pink1_sirt1.py       # 3D Atomic Contact Analyzer for PINK1 (6EQI) & SIRT1 (4ZZJ)
│   ├── analyze_2nsx.py              # 3D Atomic Contact Analyzer for GBA1 + Chaperone (2NSX)
│   └── analyze_details.py           # Detailed Sub-Pocket Geometry Inspector
├── examples/
│   ├── run_sample_patient_walkthrough.py  # Self-Contained End-to-End Sample Patient Runner (SAMPLE-PD-001)
│   └── output/                            # Generated Sample Patient ODE, CIF, & Lab Gate Telemetry
├── structures/
│   ├── 6FCX/pdb_00006fcx.cif        # Human MTHFR + FAD 701 + SAH 702 (2.50 Å)
│   ├── 1B68/pdb_00001b68.cif        # Human APOE4 Arg112-Glu109 2.53 Å Salt Bridge (2.00 Å)
│   ├── 7CK6/pdb_00007ck6.cif        # Human TOMM40 Import Pore + PC1 401 + TOMM22/TOMM7 (3.40 Å)
│   ├── 6EQI/pdb_00006eqi.cif        # Active Phosphorylated PINK1 + Ubiquitin (3.10 Å)
│   ├── 4ZZJ/pdb_00004zzj.cif        # Human SIRT1 + Carba-NAD+ + STAC 4TQ + ALY Substrate (2.74 Å)
│   └── 2NSX/pdb_00002nsx.cif        # Human GBA1 + Pharmacological Chaperone IFM 506 (2.20 Å)
├── data/                            # Empirical Ensembl VEP (GRCh38), AlphaMissense, Mutfunc, gnomAD v4, & OpenTargets JSONs
└── skills/
    └── bioenergetic-precision-protocol/
        ├── SKILL.md                 # Reusable Agent Skill Definition
        ├── scripts/bcr_pipeline.py  # Self-Contained CLI (simulate, inspect-cif, verify-gate)
        └── references/erp_ode_reference.md
```

---

## Quick Start Commands

```bash
# 1. Run the End-to-End 4-Step Sample Patient Walkthrough (SAMPLE-PD-001)
python3 examples/run_sample_patient_walkthrough.py

# 2. Run the Full 5-Cohort Biophysical ODE Simulation Benchmark
python3 simulator/sim_bioenergetic_reset.py

# 3. Use the Standalone Skill CLI (simulate / inspect-cif / verify-gate)
python3 skills/bioenergetic-precision-protocol/scripts/bcr_pipeline.py simulate \
  --cohort all --pink1-class Class_I_ATP_Pocket --days 90 --output /tmp/bcr_sim.json

python3 skills/bioenergetic-precision-protocol/scripts/bcr_pipeline.py verify-gate \
  --gdf15-day0 1950 --gdf15-current 510 --day 38 --lean-mass-delta-pct -0.6 --bhb-mmol 1.8 --output /tmp/bcr_gate.json
```

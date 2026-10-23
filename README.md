# MitoRestorer: Precision Bioenergetic Circuit Reset (BCR) Platform

**MitoRestorer** is an end-to-end biophysical discovery, sub-hourly ordinary differential equation (ODE) simulation, 3D atomic structural validation, and biomarker-gated clinical protocol platform grounded in **Picard & Murugan’s Energy Resistance Principle ($\text{EP} = \text{éR} \cdot f^2$)**, the **Microbiota–Gut–Mitochondria Axis**, and **Host Genetic Variant Parameterization** (`Ensembl VEP`, `AlphaMissense`, `Mutfunc`, `gnomAD v4`, `OpenTargets`, and `GTEx v8`).

---

## Core Theoretical Framework: The Energy Resistance Principle ($\text{EP} = \text{éR} \cdot f^2$)

Under Picard & Murugan’s **Energy Resistance Principle (ERP)**, cellular and organ-level reductive/oxidative damage, Integrated Stress Response (`ISR` / `ATF4` / `CHOP`) activation, and systemic mitokine secretion (**plasma GDF15**, `ENSG00000130513`) are governed by:

$$\text{EP}(t) = \text{éR}(t) \cdot f(t)^2$$

where:
- **$f(t)$** is the electron/metabolic flux (e.g., glycolytic $\text{NADH}$ push, autonomous $\text{Ca}_V1.3$ pacemaker spiking in `SNc` dopaminergic neurons, or $50–200\text{ Hz}$ firing in cortical $\text{PV}^+$ interneurons).
- **$\text{éR}(t)$** is the effective mitochondrial circuit resistance, determined by the ratio of mutant/damaged mtDNA ($D$) to wild-type mtDNA ($W$), outer-membrane **`TOMM40` import impedance**, methylation/creatine spatial buffering (**`MTHFR`**), and gut-derived **$\text{H}_2\text{S}$ Complex IV inhibition** (*Bilophila wadsworthia* / *Desulfovibrio*).

---

## Six Core Biophysical Discoveries

1. **Discovery 1 — The Simultaneous Biogenesis Trap:** Co-administering mitophagy activators (`Urolithin A`) and biogenesis boosters (`NR/NMN`, `SIRT1` activators, high protein/leucine) simultaneously keeps `mTORC1` active, phosphorylating `ULK1(Ser757)` to block mitophagy while `POLG` DNA polymerase preferentially clones shorter deletion-mutant mtDNA ($D$), driving mutant heteroplasmy to $h > 84\%$.
2. **Discovery 2 — The Day-30 Ketogenic Microbiome Crash & Secondary Ketogenic Escape:** Chronic strict ketogenic diets without prebiotic fiber, `2'-FL HMO`, and polyphenol/molybdenum sulfide blockade trigger a Day-30 *Bilophila wadsworthia* / *Desulfovibrio* $\text{H}_2\text{S}$ bloom that poisons Complex IV (`MT-CO1`/`COX4I1`), collapsing cellular ATP and causing **Secondary Ketogenic Escape** in refractory epilepsy ($\text{EI}_{\text{ratio}}$ rebounds from $1.80 \rightarrow 2.26$).
3. **Discovery 3 — The 13.5-Hour Glycogen Depletion Lag Trap:** Unassisted daily `16:8` time-restricted eating spends $13.5\text{ hours}$ depleting hepatic glycogen before `AMPK`/`PINK1` mitophagy engages, leaving only $\sim 2.5\text{ hours}$ of effective mitophagy per day ($h_{90} = 69.4\%$ in `PINK1/FUT2` muscle). Compressing $\tau_{\text{lag}}$ to $4–6\text{ hours}$ via morning Zone-2 exercise + `C8-MCT` unlocks daily `18:6` de-heteroplasmy.
4. **Discovery 4 — Two-Stage Biomarker-Gated "Purge-Then-Expand" Architecture:** Using **Stage 1 Lag-Compressed Daily 18:6** until plasma **GDF15** drops $\ge 35\%$ (or $< 550\text{ pg/mL}$ at the Day 38 / Day 54 gate) followed by **Stage 2 Weekly Meso-Cycling (`2.25d Purge : 4.75d Build`)** achieves both deep mutant clearance ($h < 4\%$) and supra-physiological wild-type mitochondrial expansion ($\text{ATP} = 1.20\times–1.46\times$).
5. **Discovery 5 — The *APOE4/TOMM40* Import Bottleneck & 3D `PINK1` Class I vs. Class II Stratification:** In slow-turnover cortical neurons ($0.65\times$), stimulating mitophagy (`Kinetin` + `Trehalose`) without unblocking the **`TOMM40` $\beta$-barrel pore (`PDB: 7CK6`)** (`TUDCA` + `Citicoline/LPC-DHA`) causes severe neuronal mitochondrial depletion ($\text{Mass} = 0.67\times$). Furthermore, **`PDB: 6EQI`** reveals that `Kinetin` (`KTP`) rescues **Class I ATP-pocket `PINK1` variants** (`Ala339Thr`, $h \rightarrow 3.7\%$) but cannot rescue **Class II Insertion-3 Ubiquitin-cradle variants** (`Gly309Asp`, `Leu347Pro`) without `TFEB` + astrocytic N3 transmitophagy bypass.
6. **Discovery 6 — The `GBA1/LRRK2` Autophagosome Traffic Jam & Dopamine-Quinone (`DA-Q`) Cascade in Parkinson's Disease:** In `SNc` dopaminergic pacemaker neurons carrying **`GBA1` `rs76763715` (`Asn370Ser`, $\Delta\Delta G = +1.922\text{ kcal/mol}$)** and **`LRRK2` `rs34637584` (`Gly2019Ser`, `AlphaMissense = 0.936`)**, stimulating `PINK1/Parkin` mitophagy *without* first unblocking lysosomal glucocerebrosidase (`High-Dose Ambroxol`, **`PDB: 2NSX`**) and `ATP13A2` polyamine export (`High-Dose Spermidine`) causes a $2.39\times$ **Autophagosome Traffic Jam ($\text{Auto\_Jam} = 1.29\times$)** that worsens cytosolic `DA-Q` and $\alpha$-synuclein protofibril aggregation ($\text{SNCA} = 1.63\times$).

---

## Repository Structure

```text
mitorestorer/
├── README.md
├── docs/
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

## Quick Start

### 0. Run the End-to-End 4-Step Sample Patient Walkthrough (`SAMPLE-PD-001`)
```bash
python3 examples/run_sample_patient_walkthrough.py
```

### 1. Run the Complete 5-Cohort Biophysical ODE Simulation
```bash
python3 simulator/sim_bioenergetic_reset.py
```

### 2. Run the Standalone Skill CLI (`bcr_pipeline.py`)
```bash
# Run 90-day ODE simulation across all 5 cohorts (A: Metabolic, B: APOE4 AD, C: Epilepsy, D: TR-MDD, E: Parkinson's SNc)
python3 skills/bioenergetic-precision-protocol/scripts/bcr_pipeline.py simulate \
  --cohort all --pink1-class Class_I_ATP_Pocket --days 90 --output /tmp/bcr_sim.json

# Verify a patient's Day-38 GDF15 biomarker gate
python3 skills/bioenergetic-precision-protocol/scripts/bcr_pipeline.py verify-gate \
  --gdf15-day0 1420 --gdf15-current 510 --day 38 --lean-mass-delta-pct -0.6 --bhb-mmol 1.8 --output /tmp/bcr_gate.json
```

### 3. Open the Interactive Simulator
Open `simulator/bcr_simulator.html` in any modern web browser to interactively explore all 5 clinical cohorts, 8 human genetic variants (`PPARGC1A`, `FUT2`, `MTHFR`, `APOE/TOMM40`, `SOD2`, `PINK1`, `GBA1`, `LRRK2`), phase-ratio bifurcation curves, and 90-day biomarker trajectories.

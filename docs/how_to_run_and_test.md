# How to Run & Test MitoRestorer v2.0 (Local Software Verification & Real-World Wet-Lab / Apiary Testing)

This guide covers both:
1. **Part 1 — Local Software & 3D Structural Verification:** How to run the 16-point automated test suite, the 3 end-to-end sample case scripts, the interactive HTML5 simulator, and the 3D crystal structure viewer.
2. **Part 2 — Real-World Biological & Field Testing:** Step-by-step instructions and pass/fail endpoints for testing the **14-Day Honeybee (`Apis-BCR`) Hive Patty Protocol** and the **10-Day (`240-Hour`) *Ex Vivo* CAR-T (`Bioreactor-BCR`) Manufacturing Protocol**.

---

## Part 1: How to Run & Verify the Software Locally (4 Ways)

All repository files are located in `/Users/o2satz/.gemini/antigravity/scratch/mitorestorer`.

### 1. Run the 16-Point Automated Verification & Regression Test Suite (`16 PASSED, 0 FAILED`)
[`tests/run_all_verifications.py`](../tests/run_all_verifications.py) automatically verifies all **8 local 3D `.cif` crystal structures** (`6FCX`, `1B68`, `7CK6`, `6EQI`, `4ZZJ`, `2NSX`, `2Q8G`, `1SQB`) and asserts all **8 biophysical discoveries** across Cohorts A–G:
```bash
cd /Users/o2satz/.gemini/antigravity/scratch/mitorestorer
python3 tests/run_all_verifications.py
```

### 2. Run the 3 End-to-End Case Walkthrough Scripts (Terminal + JSON Outputs)
Each script executes a complete 4-step workflow (**Step 1:** Baseline Profile & 3D `.cif` Atomic Contact Check $\to$ **Step 2:** Why Standard Care Fails $\to$ **Step 3:** Phased BCR Simulation $\to$ **Step 4:** Biomarker Gate Verification & Recipe Output):
```bash
# A. Honeybee Colony Rescue (Colony APIS-CA-ALMOND-01 | 14 Days | PDB: 1SQB)
python3 examples/run_honeybee_apis_walkthrough.py

# B. Ex Vivo CAR-T Bioreactor (Batch CART-DLBCL-042 | 240 Hours | PDB: 2Q8G)
python3 examples/run_cart_bioreactor_walkthrough.py

# C. Human Early Parkinson's Disease (Patient SAMPLE-PD-001 | 90 Days | PDB: 2NSX)
python3 examples/run_sample_patient_walkthrough.py
```
* **Generated JSON Telemetry Outputs (`examples/output/`):**
  * **Honeybee (`Apis-BCR`):** [`examples/output/apis_colony_summary.json`](../examples/output/apis_colony_summary.json) & [`examples/output/apis_step1_1sqb_azo_contacts.json`](../examples/output/apis_step1_1sqb_azo_contacts.json)
  * **CAR-T (`Bioreactor-BCR`):** [`examples/output/cart_batch_summary.json`](../examples/output/cart_batch_summary.json) & [`examples/output/cart_step1_2q8g_azx_contacts.json`](../examples/output/cart_step1_2q8g_azx_contacts.json)
  * **Human Parkinson's (`PD-BCR`):** [`examples/output/sample_patient_summary.json`](../examples/output/sample_patient_summary.json)

### 3. Use the Standalone Skill CLI (`bcr_pipeline.py`)
Run custom simulations, inspect any `.cif` crystal structure, or evaluate a patient's Day-38/54 `GDF15` blood test gate directly from the command line:
```bash
# Simulate all 7 cross-species cohorts at once
python3 skills/bioenergetic-precision-protocol/scripts/bcr_pipeline.py simulate \
  --cohort all --pink1-class Class_I_ATP_Pocket --days 90 --output /tmp/bcr_all_7.json

# Inspect atomic contacts in any 3D crystal structure (.cif)
python3 skills/bioenergetic-precision-protocol/scripts/bcr_pipeline.py inspect-cif \
  --cif structures/1SQB/pdb_00001sqb.cif --target-chain C --target-res AZO --max-dist 3.8 --output /tmp/1sqb_contacts.json

# Verify the Day-38 Plasma GDF15 Lab Gate
python3 skills/bioenergetic-precision-protocol/scripts/bcr_pipeline.py verify-gate \
  --gdf15-day0 1950 --gdf15-current 510 --day 38 --lean-mass-delta-pct -0.6 --bhb-mmol 1.8 --output /tmp/gate_check.json
```

### 4. Launch the Interactive Web Simulator & 3D Structure Viewer
* **Web Simulator:** Open [`simulator/bcr_simulator.html`](../simulator/bcr_simulator.html) (or [`index.html`](../index.html)) in any browser:
  1. Click **`▶ Load & Run Sample Patient (Step 1)`** at the top to step through the 4-stage interactive tutorial (`Step 1 → Step 2 → Step 3 → Step 4`).
  2. Use the **Clinical Genotype Archetype** dropdown (top-right) to switch between **Cohorts A–E**, **Cohort F (`Ex Vivo CAR-T Bioreactor`)**, and **Cohort G (`Apis-BCR Honeybee Rescue`)**.
  3. Click **`Action Plan & Lab Gate`** to test custom blood labs in the interactive Traffic-Light calculator, or click **`Disease Circuit & 3D Atlas`** to inspect all 7 cohorts and 8 crystal structures side-by-side.

---

## Part 2: How to Run & Test for Real-World Biological Results

### Track 1: How to Run the 14-Day Honeybee (`Apis-BCR`) Apiary or Cage Test
*Full formulation manual:* [`docs/apis_bcr_beekeeper_protocol.md`](apis_bcr_beekeeper_protocol.md)

#### 1. Experimental Setup (`3 Groups of Hives or Lab Cages` post-spray / during pollination)
* **Group 1 (Control — Standard Practice):** Feed 50% High-Fructose Corn Syrup (HFCS) + standard pollen substitute for 14 days.
* **Group 2 (Unphased Vitamin Control):** Feed Stage 2 CoQ10/Vitamin patty continuously from Day 1 without Stage 1 detox/mitophagy.
* **Group 3 (`Apis-BCR` 2-Stage Test Group):**
  * **Days 1–5 (`Patty A — In-Hive Detox & Mitophagy Patty`, `250 g/hive`):** Sucrose/Trehalose base (*zero HFCS*) + `120 mg/kg` ***p*-coumaric acid** (or `8 g/kg` propolis powder; upregulates `CYP9Q1/2/3` to clear strobilurins from `PDB: 1SQB`) + `150 mg/kg` **spermidine** + `250 mg/kg` **Urolithin A** + `2 g/kg` live ***Snodgrassella alvi***.
  * **Days 5–14 (`Patty B — Forager OXPHOS Build Patty`, swap at Day 5, `350 g/hive`):** `400 mg/kg` **water-dispersible Ubiquinol (`CoQ10`)** + `300 mg/kg` **Nicotinamide Riboside (`NR`)** + `80 mg/kg` **Riboflavin (`Vit B2`)** + `500 mg/kg` **24-methylenecholesterol**.

#### 2. Pass/Fail Endpoints to Measure (Day 5 & Day 14)

| Assay / Measurement | Method | Expected Control Outcome (Group 1 HFCS) | Pass Target (`Apis-BCR` Group 3) |
| :--- | :--- | :---: | :---: |
| **1. Forager Homing Survival (Day 10–14)** | Paint-mark or RFID-tag `100 foragers/hive`, release `1.5 km` away, count 24h return | $18.0\text{–}30.0\%$ | **$\ge 90.0\%$ (Model: $96.5\%$)** |
| **2. 60-Min Flight Endurance** | Tethered flight-mill distance ($\text{km}$) & optical wingbeat tachometer ($\text{Hz}$) | Rapid fatigue ($<15\text{ min}$) | **Sustained $100\text{–}250\text{ Hz}$ flight** |
| **3. Thoracic Flight-Muscle ATP** | Bioluminescent luciferase ATP/ADP kit on 5 homogenized thoraxes | $0.06\text{–}0.25\times$ | **$\ge 1.50\times$ (Model: $1.74\times$)** |
| **4. Complex III $Q_o$ & *S. alvi* Recovery** | Complex III oxidoreductase assay + hindgut *S. alvi* 16S qPCR | $72.9\%$ $Q_o$ block, $0.16\times$ *S. alvi* | **$\le 10\%$ $Q_o$ block, $\ge 0.90\times$ *S. alvi*** |

---

### Track 2: How to Run the 10-Day (`240-Hour`) CAR-T Bioreactor Bench Test
*Full GMP media recipe & sensor SOP:* [`docs/cart_bioreactor_sop.md`](cart_bioreactor_sop.md)

#### 1. Experimental Setup (`3 Parallel Wells in a 24-Well G-Rex Plate` or Perfusion Bioreactor)
* **Well 1 (Standard Control):** `25 mM` Glucose + `IL-2` (`100 IU/mL`) for 240 hours.
* **Well 2 (Static Mitophagy Control):** `Urolithin A` (`25 μM`) + `Rapamycin` (`10 nM`) continuously for 240 hours.
* **Well 3 (`Bioreactor-BCR` Biphasic Test):**
  * **Hours 0–64 (`Media-A — Purge & Unmasking`):** `4 mM` Glucose + `10 mM` Galactose + `2 mM` D-BHB + **`1 μM` `AZD7545` (or `2 mM` `DCA`, blocks `PDK1` `PDB: 2Q8G`)** + **`25 μM` Urolithin A** + **`10 μM` Kinetin (`PDB: 6EQI`)** + `50 μM` Citicoline + `IL-7/IL-15` (`10 ng/mL`).
  * **Hour 64 Gate Check:** Verify supernatant `Lactate:Pyruvate` $\le 13.0$ and `GDF15` $< 380\text{ pg/mL}$, then wash and switch to **`Media-B` (Hours 64–240 — Clonal Expansion):** `12 mM` Glucose + **`500 μM` NR** + **`5 μM` Resveratrol (`PDB: 4ZZJ`)** + **`100 μM` Citicoline (`PDB: 7CK6`)** + **`2 mM` Creatine** + **`4 mM` L-Arginine** + `IL-7/IL-15` (`15 ng/mL`).

#### 2. Pass/Fail Release Endpoints at Hour 240 (Day 10 Harvest)

| Assay / Release Metric (`Hour 240`) | Well 1 (`High-Glc + IL-2`) | Well 2 (`Static UroA + Rapa`) | Pass Target (`Bioreactor-BCR` Well 3) |
| :--- | :---: | :---: | :---: |
| **1. Total Clonal Fold Expansion** | $51.9\times$ | $10.8\times$ *(FAIL $<45\times$)* | **$\ge 45.0\times$ (Model: $62.9\times$ PASS)** |
| **2. $\text{T}_{\text{SCM}}$ Stemness (`CD45RA+ CCR7+ CD62L+ TCF7+`)** | $10.4\%$ *(FAIL)* | $72.8\%$ | **$\ge 60.0\%$ (Model: $78.0\%$ PASS)** |
| **3. Depolarized Mitochondria (`TMRE-low / MitoTracker+`)** | $56.1\%$ *(FAIL)* | $0.2\%$ | **$\le 10.0\%$ (Model: $1.6\%$ PASS)** |
| **4. Exhaustion (`PD-1+ TIM-3+ TOX+`)** | $60.2\%$ *(FAIL)* | $4.9\%$ | **$\le 15.0\%$ (Model: $4.8\%$ PASS)** |
| **5. Seahorse Spare Respiratory Capacity (`SRC`)** | $0.93\times$ *(FAIL)* | $2.69\times$ | **$\ge 2.50\times$ (Model: $3.78\times$ PASS)** |

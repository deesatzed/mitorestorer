# MitoRestorer Step-by-Step User Walkthrough (Executive & Clinical Mechanistic Guide)

Welcome to **MitoRestorer**! This guide walks you step-by-step from zero background knowledge to running a complete 90-day cellular simulation, understanding why standard supplement stacks fail, and using the **Day-38 Blood Test Gate** to generate a personalized daily protocol.

---

## The Big Picture in 60 Seconds (The "Rechargeable Battery" Analogy)

Every cell in your body—especially high-energy organs like your **brain**, **heart**, and **skeletal muscles**—is powered by thousands of tiny rechargeable batteries called **mitochondria**.

Over time, or because of inherited genetic variants (SNPs), two things go wrong:
1. **Broken Batteries Pile Up ($h\%$ / Heteroplasmy):** Some mitochondria suffer DNA damage ($D$). Instead of producing clean energy (`ATP`), they leak toxic free radicals and raise **Energy Resistance ($\text{éR}$)** inside the cell.
2. **Pushing a Clogged Circuit Causes Cellular Burnout ($\text{EP} = \text{éR} \cdot f^2$):** Under Picard & Murugan’s **Energy Resistance Principle**, when you force high metabolic flux ($f$, such as high sugar intake or continuous neural firing) through a high-resistance mitochondrial circuit ($\text{éR}$), cellular stress ($\text{EP}$) increases with the **square** of the flux ($f^2$). The cell sounds a systemic "smoke alarm" in the blood called **GDF15**.

**The Goal of MitoRestorer:**
Instead of pushing a clogged circuit or taking random supplements all at once, MitoRestorer calculates how to **clean out the broken batteries first (Stage 1 Purge)** and **copy the healthy batteries second (Stage 2 Build)**—guided by your specific DNA variants, gut microbiome, and 3D protein structures.

---

## Step 1 of 4: Choose a Patient Scenario & Inspect Their Genetic Bottlenecks

### How to Do It:
* **In the Web App (`simulator/bcr_simulator.html`):** Look at the top-right dropdown labeled **"Clinical Genotype Archetype"** (or click the green **`▶ Load & Run Sample Patient (Step 1)`** button in the top Guided Walkthrough bar).
* **In the Terminal:** Run `python3 examples/run_sample_patient_walkthrough.py`.

### What You Are Choosing:
You can pick from **5 pre-calibrated clinical scenarios** (or build a custom combination using the **Empirical SNP Panel** on the left):
1. **Cohort A (Metabolic Syndrome / Muscle Fatigue):** Carrying `PINK1` (slow battery trash collector) + `FUT2 rs601338 A/A Non-Secretor` (a gut microbiome variant present in **22.5% of people** that prevents gut bacteria from turning pomegranate/berries into the battery-cleaning postbiotic **Urolithin A**).
2. **Cohort B (Alzheimer’s / Cognitive Risk — `APOE4 / TOMM40`):** Carrying `APOE4` (`rs429358`), which clogs the **`TOMM40` import door** (`PDB: 7CK6`) on the outside of neuronal mitochondria.
3. **Cohort C (Refractory Epilepsy — `Epi-BCR`):** Fast-firing inhibitory brain cells ($\text{PV}^+$ interneurons) that suffer "Secondary Ketogenic Escape" after Day 30 on a classic ketogenic diet when gut bacteria (*Bilophila* / *Desulfovibrio*) produce **hydrogen sulfide ($\text{H}_2\text{S}$)** gas that poisons Complex IV.
4. **Cohort D (Treatment-Resistant Depression — `Affective-BCR`):** Gut-driven microglial inflammation (**Quinolinic Acid / `QUIN`**) and high **GDF15** lock the brain into an energy-conservation shutdown, causing rapid relapse after standalone ketamine infusions if mitochondrial resistance ($\text{éR}$) isn't cleared first.
5. **Cohort E (Parkinson’s Disease — `PD-BCR` / Sample Case `SAMPLE-PD-001`):** Dopamine-producing pacemaker neurons in the *substantia nigra* (`SNc`) carrying **`GBA1 rs76763715`** (`Asn370Ser`) and **`LRRK2 rs34637584`** (`Gly2019Ser`), where the cell's lysosomal "incinerator" is jammed (`58%` damaged batteries at Day 0).

---

## Step 2 of 4: Test Why Standard Treatments & "All-at-Once" Supplements Backfire

### How to Do It:
* Click **`Next Step →`** in the top Guided Walkthrough bar (or click **Regimen 4: Unphased "All-in-One" Stack** or **Regimen 3: Chronic Strict Keto Alone** in Box 1 on the left).

### What to Look For on the Charts:
Look at **Chart A (Mitochondrial Selection)** and the top-left **Damaged mtDNA ($h\%$)** KPI card:
* **Why Regimen 4 (Taking Everything Together) Fails — *Discovery 1 (The Simultaneous Biogenesis Trap)*:**
  * If a patient takes "cleaning" supplements (`Urolithin A`) and "building" supplements (`NAD+ / NR` + high protein) **at the same time every day**, the cell's growth sensor (`mTORC1`) stays turned ON.
  * `mTORC1` locks the trash collector (`ULK1`) shut while telling the DNA copying enzyme (`POLG`) to copy existing mitochondria.
  * Because broken mitochondrial DNA loops ($D$) have missing sections, they are physically shorter and **copy 32% faster** than healthy wild-type loops ($W$). Taking everything together accidentally **clones the broken batteries to >84%**!
* **Why Pushing Mitophagy Alone in `GBA1` Parkinson’s Backfires — *Discovery 6 (The Autophagosome Traffic Jam)*:**
  * In **Cohort E (Parkinson’s)**, if you trigger battery-cleaning (`Urolithin A` + `Kinetin`) while the cell's lysosomal incinerator (`GBA1`) is trapped in the endoplasmic reticulum (`PDB: 2NSX`), wrapped-up trash bags (**mitophagosomes**) pile up inside the axon (`Auto_Jam = 1.29x`), leaking toxic **Dopamine-Quinone (`DA-Q`)** and worsening $\alpha$-synuclein clumping!

---

## Step 3 of 4: Activate the Phased Bioenergetic Circuit Reset (BCR) & 3D Cofactor Rescue

### How to Do It:
* Click **`Next Step →`** (to **Step 3 of 4**) or select **Regimen 5: Biphasic Chrono-Synbiotic (BCR)** on the left.
* Click the **`Disease Circuit & 3D Atlas`** tab to see how each 3D protein crystal structure guides the exact supplement rescue.

### How the 3-Stage Reset Solves the Problem:
1. **Stage 0 (Days 1–14) — Unblock the Doors & Incinerators First:**
   * **For `APOE4/TOMM40` (`PDB: 7CK6`):** Use **TUDCA + Citicoline/LPC-DHA** to reopen the mitochondrial protein import pore.
   * **For `GBA1/LRRK2` Parkinson’s (`PDB: 2NSX`):** Use **High-Dose Ambroxol + Spermidine (`ATP13A2` substrate)** to chaperone `GBA1` into the lysosome *before* starting mitophagy.
   * **For the Gut Microbiome (`FUT2` & $\text{H}_2\text{S}$):** Use **Tributyrin + 2'-FL HMO + Molybdenum + Polyphenols** to block toxic $\text{H}_2\text{S}$ gas.
2. **Stage 1 (Days 15–38) — Lag-Compressed Daily 18:6 Cleaning Purge:**
   * Why does standard 16:8 intermittent fasting often fail? Because the liver takes **13.5 hours** to burn through stored glycogen sugar before cellular cleaning (`PINK1` mitophagy) even starts—leaving only 2.5 hours of cleaning!
   * By combining a **30-minute morning Zone-2 walk + C8-MCT Oil** with **Direct Urolithin A (`1,000 mg`) + Kinetin (`PDB: 6EQI`)**, MitoRestorer compresses that lag from `13.5 hours` down to **`5–6 hours`**, giving the cell **12+ hours of deep cleaning every day** while holding growth boosters (`NR/NMN`) OFF.
3. **Stage 2 (Days 38–90) — Weekly Anabolic Rebuilding (`2.25d Purge : 4.75d Build`):**
   * Once broken batteries ($D$) are cleared below `8%`, the protocol switches to **4.75 days/week of rebuilding** using **NR/NMN + *PncA* Probiotic + Resveratrol (`SIRT1 PDB: 4ZZJ`) + Riboflavin B2 (`MTHFR PDB: 6FCX`) + Creatine (`5 g`)**. Now, when the cell copies its batteries, almost 100% of the templates are healthy wild-type ($W$) mitochondria!

---

## Step 4 of 4: Verify the Day-38 Blood Test Gate & Generate the Daily Action Plan

### How to Do It:
* Click **`Next Step →`** (to **Step 4 of 4**) or open the **`Action Plan & Lab Gate`** tab in the simulator.
* Click **`Load Sample Patient Day-38 Labs`** (or enter any custom patient lab numbers).

### How the Traffic-Light Gate Works:
In real life, how do you know when a patient is ready to switch from **Stage 1 (Cleaning)** to **Stage 2 (Rebuilding)**? You check 3 simple numbers at **Day 38**:
1. **Fasting Plasma GDF15 ("The Cellular Smoke Alarm"):**
   * **🟢 GREEN LIGHT:** Dropped by **$\ge 35\%$** from Day 0 (or is **$< 550\text{ pg/mL}$**). In our sample patient `SAMPLE-PD-001`, GDF15 dropped from `1,950 pg/mL` to `510 pg/mL` (**`-73.8%`**), triggering an immediate Green Light to start Stage 2!
   * **🟡 YELLOW LIGHT:** Dropped by `< 35%` at Day 38 with stable muscle mass $\rightarrow$ stay in Stage 1 for 16 more days until Day 54.
2. **Morning Capillary Ketones ($\beta\text{-OHB}$ at 09:00 AM):**
   * Target **`1.2 – 2.8 mmol/L`**, confirming that the morning walk + C8-MCT successfully compressed the glycogen lag.
3. **Muscle Mass Safety Override (DEXA / BIA):**
   * **🔴 RED OVERRIDE:** If a patient loses **$\ge 2.5\%$ lean muscle mass**, they immediately switch to Stage 2 Anabolic Build regardless of GDF15.

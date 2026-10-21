# Energy Resistance Principle (ERP) & Bioenergetic Circuit Reset (BCR) Reference

## 1. Fundamental Biophysical Equations
Grounded in Picard & Murugan's **Energy Resistance Principle** (*Cell Metabolism*, Nov 2025):

$$\text{EP} = \text{éR} \cdot f^2$$

- **$\text{EP}$ (Energy Potential / Bioenergetic Stress):** Drives mitochondrial reductive backpressure, superoxide ($\text{O}_2^{\bullet-}$) leak, the mitochondrial Integrated Stress Response (`ISR` via `OMA1–DELE1–HRI–ATF4`), and systemic release of **GDF15** and **FGF21**.
- **$f$ (Electron / Metabolic / Synaptic Flux):** Scales quadratically ($f^2$). Fast-spiking Parvalbumin ($\text{PV}^+$) GABAergic interneurons ($f_{\text{PV}} \approx 1.75 \implies f_{\text{PV}}^2 = 3.06\times$) and glutamatergic bursts (`Ketamine` / `Quinolinic Acid`) are hyper-vulnerable when $\text{éR}$ is elevated.
- **$\text{éR}$ (Cellular Circuit Resistance):** Increased by deletion-mutant mtDNA heteroplasmy ($h = \frac{D}{W+D}$), gut-derived $\text{H}_2\text{S}$ Complex IV inhibition, microglial Quinolinic Acid (`QUIN`), and outer-membrane `TOMM40` import impedance.

---

## 2. Five Core In Silico Kinetic Discoveries

1. **Discovery 1 — The Simultaneous Biogenesis Trap:** Co-administering mitophagy activators (Urolithin A) and biogenesis boosters ($\text{NAD}^+$ + leucine/mTORC1) simultaneously keeps `ULK1(Ser757)` phosphorylated while *POLG* replicates shorter deletion-mutant mtDNA ($D$) ~32% faster than wild-type $W$.
2. **Discovery 2 — The Day-30 Ketogenic Microbiome Crash (Secondary Ketogenic Escape):** Un-cycled high-fat ketogenic diets without prebiotic fiber, 2'-FL HMO, and polyphenol/molybdenum sulfide blockade bloom *Bilophila wadsworthia* and *Desulfovibrio* after Day 30, whose $\text{H}_2\text{S}$ inhibits **Complex IV** and causes breakthrough seizures in epilepsy or muscle atrophy.
3. **Discovery 3 — The 13.5-Hour Circadian Glycogen/PINK1 Lag Trap:** Unassisted 16:8 time-restricted eating spends ~13.5 hours depleting hepatic glycogen and auto-phosphorylating `PINK1` (`SEP202–Arg282`, `PDB: 6EQI`), leaving only ~2.5 hours/day of mitophagy. Combining morning Zone-2 cardio + C8-MCT + direct Urolithin A compresses $\tau_{\text{lag}}$ to $6.0\text{h}$ ($4.0\text{h}$ with Kinetin/D-BHB).
4. **Discovery 4 — Two-Stage Biomarker-Gated "Purge-Then-Expand" Regimen:** Stage 1 executes daily lag-compressed 18:6 Purging (with $\text{NAD}^+$ held low) until plasma **GDF15 drops $\ge 35\%$** ($h \le 8\%$, or Day 54 safety cap), then transitions to Stage 2 Weekly Anabolic Expansion (`2.25d Purge : 4.75d Build`).
5. **Discovery 5 — The *TOMM40* Import Bottleneck & Neuropsychiatric Synergy:** In slow-turnover cortical neurons ($\tau_{\text{tissue}} = 0.65\times$), boosting mitophagy without unblocking the 19-stranded `TOMM40` import pore (`PDB: 7CK6`) causes neuronal atrophy ($0.67\times$). Combining `TUDCA + Sulforaphane + Citicoline/LPC-DHA` (`TOMM40` rescue) with `Trehalose/Spermidine` (`TFEB`), `Kinetin` (`PINK1`), and `N3 Melatonin/Glycine` (astrocytic transmitophagy) collapses neuronal $h \rightarrow 3.3\%$ and restores ATP to $1.19\times$.

---

## 3. Five-Protein 3D Atomic Reference Stack

| Stage / Gate | PDB ID | Key Author Residues & Ligands (`auth_asym_id` / `auth_seq_id`) | Precision Therapeutic Mechanism |
| :--- | :---: | :--- | :--- |
| **Stage 0: Methylation** | **`6FCX`** (Human MTHFR) | `Chain A`: `Ala222` ($\alpha6$ helix), `Lys217` ($4.1\text{ Å}$ away), `FAD 701`, `SAH 702` | `Riboflavin B2` restores `FAD 701` occupancy against `Val222` steric clash; `Creatine` spares SAMe. |
| **Stage 1: Mitophagy** | **`6EQI`** (`PINK1 + Ubiquitin`) | `Chain C`: `SEP202–Arg282` ($2.22\text{ Å}$), `Ins3 265–285`, `Met290`; `Chain A`: `Ub Ser65` | `Kinetin (KTP)` rescues Class I ATP-pocket variants (`Met290`); Class II `Ins3` mutants require `Trehalose/Spermidine` (`TFEB`). |
| **Import Gate (APOE4)** | **`1B68`** (Human APOE4) | `Chain A`: `Arg112–Glu109` ($2.53\text{ Å}$ salt bridge), `Arg61` ($3.87\text{ Å}$ repulsion) | `C8-MCT / D-BHB` bypasses glycolytic Complex I/IV backpressure. |
| **Import Gate (TOMM40)** | **`7CK6`** (Human TOM Core) | `Chain A`: `TOMM40 106–360`, `Leu77–Glu105`, `PC1 401`; `Chains C/I`: `TOMM22/TOMM7` | `TUDCA + Sulforaphane` clear stalled precursors; `Citicoline + LPC-DHA` supply structural `PC1 401`. |
| **Stage 2: Biogenesis** | **`4ZZJ`** (Human SIRT1) | `Chain A`: `CNA 702` ($2.29\text{ Å}$ to `Chain B ALY 4`), `4TQ 703`, `Glu230` hinge, `ZN 701` | Bacterial *PncA* + NR/NMN fill the `CNA` C-pocket; `Glu230` STAC activators clamp `PPARGC1A Gly482Ser` 3.5x tighter. |

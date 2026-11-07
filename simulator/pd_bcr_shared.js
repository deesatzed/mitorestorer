/**
 * MitoRestorer v3.0 — Shared Parkinson's Disease (PD-BCR) 10-State ODE Engine & Suite Navigation
 * ==============================================================================================
 * Powers the 4-Page Academic PD-BCR Suite:
 *   Page 1: parkinsons_simulator.html     — Pathophysiology, Selective Mitophagy & 4-Phase Protocol
 *   Page 2: parkinsons_patients.html      — Patient Genetics, Dosing Titration & 10-State Telemetry
 *   Page 3: parkinsons_trials.html        — Retrospective Clinical Trial Graveyard Audit (10 Trials)
 *   Page 4: parkinsons_breakthroughs.html — 2025–2026 Breakthroughs, Sensitivity Sweeps & 3D Atlas
 */

(function (global) {
  'use strict';

  const STORAGE_KEY = 'mitorestorer_pd_bcr_state_v3';

  const DEFAULT_STATE = {
    genotype: 'gba1_lrrk2',
    activeRegimen: 'pd_bcr',
    gateDay: 27,
    ambroxol: 0.95,
    hypoxiaShield: 0.94,
    kinetin: 0.92,
    stage2Build: 0.96,
    vdacShield: 0.92,
    deagingEnabled: true,
    creatineMode: 'stage2_gated',
    h0: 0.58,
    genoJamScale: 1.0,
    genoHyperoxiaScale: 1.0,
    genoMitophagyScale: 1.0,
    stepIdx: 135, // Day 27.0 (0.2d per step, 450 steps = 90.0d)
    camPreset: 0,
    rotYaw: -0.62,
    rotPitch: 0.42
  };

  const GENOTYPE_PROFILES = {
    sporadic: {
      name: 'Sporadic Idiopathic PD (Complex I / Environmental)',
      badge: 'Baseline h₀ = 52.0% · Moderate Hyperoxia',
      h0: 0.52,
      gateDay: 25,
      genoJamScale: 0.85,
      genoHyperoxiaScale: 0.95,
      genoMitophagyScale: 1.0,
      snps: 'Wild-Type GBA1/LRRK2 · Somatic Complex I mtDNA deletions',
      clinicalNote: 'Standard idiopathic SNc pacemaker bioenergetic decline. Responds rapidly to Stage 0 Ambroxol + GBT601/PT2399 oxygen clamp, reaching the GDF15 < 650 pg/mL safety gate by Day 25.'
    },
    gba1_lrrk2: {
      name: 'GBA1 (N370S / rs76763715) + LRRK2 (G2019S)',
      badge: 'High Lysosomal Jam (1.15×) · Baseline h₀ = 58.0%',
      h0: 0.58,
      gateDay: 27,
      genoJamScale: 1.15,
      genoHyperoxiaScale: 1.05,
      genoMitophagyScale: 0.92,
      snps: 'GBA1 p.Asn409Ser (2NSX) + LRRK2 p.Gly2019Ser (7CK6)',
      clinicalNote: 'Severe glucocerebrosidase traffic jam and hyperactive LRRK2 kinase arrest autophagosome-lysosome fusion. Requires mandatory 14-day Stage 0 high-dose Ambroxol (95%) before initiating Stage 1 mitophagy.'
    },
    pink1_prkn: {
      name: 'PINK1 / PRKN Loss-of-Function (Early-Onset Autosomal Recessive)',
      badge: 'Impaired Mitophagy Tagging · Baseline h₀ = 62.0%',
      h0: 0.62,
      gateDay: 31,
      genoJamScale: 0.90,
      genoHyperoxiaScale: 1.10,
      genoMitophagyScale: 0.72,
      snps: 'PINK1 Kinase Domain Variant (6EQI) · PRKN Ubiquitin E3 Ligase',
      clinicalNote: 'Defective Ser65-ubiquitin phosphorylation on depolarized mitochondria slows Stage 1 clearance. Requires high Kinetin riboside (KTP, 95%) neo-substrate dosing and extending the Stage 1 purge gate to Day 31.'
    },
    apoe4_snca: {
      name: 'SNCA Triplication / High α-Synuclein + Hyperoxia Phenotype',
      badge: 'Severe SNCA & DA-o-Quinone Burn · Baseline h₀ = 60.0%',
      h0: 0.60,
      gateDay: 29,
      genoJamScale: 1.22,
      genoHyperoxiaScale: 1.20,
      genoMitophagyScale: 0.90,
      snps: 'SNCA Locus Duplication/Triplication + EPAS1 (5UFP) O₂ Mismatch',
      clinicalNote: 'High baseline α-synuclein oligomerization clogs TOMM40 import pores, driving extreme midbrain hyperoxia (pO₂ > 1.70×) and dopamine-o-quinone toxicity. Requires full GBT601 + PT2399 oxygen clamping.'
    },
    advanced_multihit: {
      name: 'Advanced Late-Stage Multi-Hit PD (High Heteroplasmy & eAge)',
      badge: 'Critical Baseline h₀ = 68.0% · Extended Day-34 Gate',
      h0: 0.68,
      gateDay: 34,
      genoJamScale: 1.25,
      genoHyperoxiaScale: 1.22,
      genoMitophagyScale: 0.82,
      snps: 'GBA1 + MTHFR C677T (6FCX) + SOD2 V16A (1B68) + PPARGC1A (4ZZJ)',
      clinicalNote: 'Near-cliff heteroplasmy (68% mutant D, only 32% healthy W). Must strictly enforce Stage 1 C8-MCT ketone ATP support and wait until Day 34 for h(t) < 8% before triggering Stage 2 PGC-1α/Creatine biogenesis.'
    }
  };

  const PHASE_DEFINITIONS = [
    {
      id: 0,
      shortTitle: 'Stage 0: Prep & Shield',
      daysLabel: 'Days 1–14',
      dayRange: [0, 14],
      headline: 'Stage 0 (Days 1–14): Unclog Lysosomal Drains & Clamp Brain Hyperoxia Before Any Purge',
      interventions: 'Ambroxol (1,200 mg/d, GBA1 2NSX) + GBT601/PT2399 (pO₂ Clamp, 5UFP) + Butyrate Synbiotic',
      survivalMechanism: 'NO mitophagy purge yet. Empties and acidifies clogged GBA1 lysosomes (Jam 0.82× → 0.15×) and clamps excess tissue O₂ (pO₂ 1.63× → 0.95×) so the neuron will not rupture or burn when mitophagy starts on Day 15.',
      purpose: 'Clear the autophagosome-lysosome bottleneck (GBA1, PDB 2NSX) and quench cytoplasmic Dopamine-o-Quinone formation by clamping midbrain hyperoxia (EPAS1/HIF-2α, PDB 5UFP).',
      whyFlags: 'Starting PINK1 mitophagy while GBA1 lysosomes are blocked causes undigested mitochondrial cargo to pile up (Jam > 1.25×), rupturing lysosomes and spiking α-synuclein.',
      decisionScope: 'ACTIVATE: Ambroxol + GBT601/PT2399 + Butyrate/TUDCA.\nHOLD OFF: Kinetin/Urolithin A mitophagy (wait until Day 15) and NAD+/Creatine biogenesis (wait until Gate).',
      methodsUsed: 'dJam/dt = k_in·u_kin·(1 - 0.8·u_amb) - (0.28 + 0.74·u_amb)·Jam\npO₂(t) → (1 + 0.95·max(0, 0.88 - OXPHOS))·(1 - 0.44·u_hyp)',
      overallImpact: 'Reduces lysosomal traffic jam below 0.20× and drops pO₂ to normoxic ~0.95× by Day 14, creating a safe runway for Stage 1 selective mitophagy.'
    },
    {
      id: 1,
      shortTitle: 'Stage 1: Selective Purge',
      daysLabel: 'Days 15–27',
      dayRange: [15, 26.8],
      headline: 'Stage 1 (Days 15–27): Voltage-Selective Mitophagy of Damaged mtDNA + Ketone ATP Survival Bridge',
      interventions: 'Kinetin Riboside (PINK1 6EQI) + Urolithin A + 18:6 Pulsed Fasting + C8-MCT Ketones (β-OHB)',
      survivalMechanism: 'Why doesn\'t purging 58% damaged mitochondria starve the neuron? (1) PINK1 only tags depolarized (low-ΔΨm) damaged D mitochondria, sparing healthy W mitochondria. (2) C8-MCT ketones boost ΔG_ATP in the remaining 42% healthy W pool so cellular ATP stays safely above the 0.45× survival floor!',
      purpose: 'Selectively cull depolarized mutant mtDNA copies (D) from 58% down to <8% while supporting neuronal pacemaker ATP above the survival floor via C8-MCT ketone oxidation.',
      whyFlags: 'Damaged D templates have a replicative advantage (r_D = 1.18×) due to deletion size and local ROS. If not purged below 8% first, any anabolic biogenesis will preferentially clone mutant D.',
      decisionScope: 'ACTIVATE: Kinetin (PINK1 Ser65-Ub activator) + Urolithin A + C8-MCT Ketone Bridge.\nHOLD OFF: High-dose NAD+/NR and Creatine (prevent premature biogenesis or Cys283 uMtCK oxidation).',
      methodsUsed: 'dD/dt = rep_D - (0.015 + 0.088·u_kin·lysoFlux)·D\nATP(t) = clip(0.86·OXPHOS + 0.22·KetoneBridge - 0.14·Jam, 0.18, 1.55) > ATP_crit (0.45×)',
      overallImpact: 'Drives mutant heteroplasmy h(t) = D/(W+D) from 58.0% down to <7.5% by Day 27 while keeping SNc ATP at ~0.78× (well above the 0.45× neuronal death floor).'
    },
    {
      id: 2,
      shortTitle: '⚡ Day-27 Safety Gate',
      daysLabel: 'Day 27 Lab Check',
      dayRange: [26.8, 28.5],
      headline: 'Biomarker Safety Gate (Day 27): Hard Brake on Mitophagy Before Over-Purging Starves the Cell',
      interventions: 'Blood Lab Verification: Plasma GDF15 < 650 pg/mL · Exosomal h(t) < 8% · α-KG:Succinate > 1.10',
      survivalMechanism: 'Continuing aggressive mitophagy after D is cleared would erode healthy mitochondrial mass and starve the pacemaker neuron. The Day-27 Biomarker Gate shuts OFF aggressive mitophagy (tapering Kinetin -70%) and flips the switch to Stage 2 rebuilding.',
      purpose: 'Verify via non-invasive blood biomarkers that mutant heteroplasmy is <8% and cellular stress (GDF15) has resolved before unlocking PGC-1α mitochondrial cloning.',
      whyFlags: 'Proposition 1 Bifurcation Boundary: Switching too early (e.g., Day 14, when h > 25%) re-expands mutant D; never switching (unbraked purge) depletes total mitochondrial mass.',
      decisionScope: 'GATE PASS (h < 8%, GDF15 < 650 pg/mL): Taper Kinetin by 70%; unlock Stage 2 NAD+ + Creatine + α-KG.\nGATE HOLD (if h ≥ 8%): Extend Stage 1 purge by 4–7 days (as in Advanced Multi-Hit PD, Day 34).',
      methodsUsed: 'Gate Condition: h(t_sw) = D(t_sw)/[W(t_sw)+D(t_sw)] ≤ h_crit (0.08) ∧ α-KG/Succinate ≥ 1.10',
      overallImpact: 'Locks in a >92% pure wild-type mitochondrial template pool at the exact inflection point for safe clonal expansion.'
    },
    {
      id: 3,
      shortTitle: 'Stage 2: Clonal Rebuild',
      daysLabel: 'Days 28–90',
      dayRange: [28.5, 90],
      headline: 'Stage 2 (Days 28–90): Clonal Wild-Type Biogenesis, uMtCK Creatine Shuttle & TET Epigenetic De-Aging',
      interventions: 'NAD⁺ Precursors + SIRT1/PGC-1α (4ZZJ) + Creatine (uMtCK 1QK1) + α-KG + 5-MTHF (MTHFR 6FCX)',
      survivalMechanism: 'Because >92% of remaining mitochondria are now healthy W templates, activating PGC-1α clones healthy mitochondria exponentially back to 108% mass, restores the octameric Creatine shuttle (97.0%), and rolls back post-mitotic epigenetic age (eAge 78% → 2.6%).',
      purpose: 'Clone the purified wild-type mitochondrial fleet (W → 1.08×), reactivate the octameric uMtCK phosphocreatine energy shuttle (PDB 1QK1), and erase aberrant DNA methylation via TET1–3 (Berg et al. 2026).',
      whyFlags: 'Now that ROS and DA-o-Quinone are quenched, uMtCK Cys283 remains unoxidized (97% octamer integrity), allowing Creatine to buffer high-frequency CaV1.3 pacemaker spikes.',
      decisionScope: 'ACTIVATE FULL ANABOLIC STACK: NR/NMN + Pterostilbene (SIRT1 4ZZJ) + Creatine Monohydrate (1QK1) + Ca-α-KG + 5-MTHF (6FCX).\nTAPER: Kinetin maintenance (30%) + Hypoxia shield taper (45%) as W respiration normalizes O₂ uptake.',
      methodsUsed: 'dW/dt = (0.024 + 0.068·u_nad·SIRT1)·W·(1 - (W+D)/1.18)\ndeAge/dt = -k_TET·(αKG/Succ)·5MTHF·eAge (Berg et al. 2026, DOI: 10.1111/acel.70751)',
      overallImpact: 'Achieves 99.1% Tyrosine Hydroxylase (TH) motor rescue, 2.4% final mutant mtDNA h(90d), 97.0% uMtCK shuttle activity, and 2.6% biological epigenetic age.'
    }
  ];

  const TRIAL_AUDIT_PRESETS = {
    steady_pd3: {
      name: '1. STEADY-PD III (Isradipine CaV1.3 Blocker)',
      nct: 'NCT02168842 · Phase III (N = 336, 36 Months)',
      citation: 'Parkinson Study Group STEADY-PD III Investigators, Ann Intern Med 2020',
      publishedResult: 'FAILED Primary Endpoint (ΔUPDRS -0.27 pts, p = 0.85)',
      bottleneck: 'Lowered pacemaker frequency f² slightly (~15%), but left baseline mutant mtDNA (h = 58%), GBA1 lysosomal jam, and brain hyperoxia untouched.',
      rescueMechanism: 'Pairing CaV1.3 pacemaker unloading with Stage 0 GBA1 unclogging and Stage 1 PINK1 mitophagy clears the underlying eR resistance.',
      gateDay: 27, ambroxol: 0.10, hypoxiaShield: 0.08, kinetin: 0.10, stage2Build: 0.15, creatineMode: 'none',
      trialCurveRegimen: 'levodopa',
      msg: 'STEADY-PD III (Isradipine, N=336) failed because partial CaV1.3 L-type calcium channel blockade only trims pacemaker demand f² by ~15% without clearing the 58% mutant mtDNA resistance (eR) or clamping midbrain hyperoxia (pO₂ = 1.61×).'
    },
    sure_pd3: {
      name: '2. SURE-PD3 (Inosine Urate Elevation)',
      nct: 'NCT02642393 · Phase III (N = 298, Halted Early)',
      citation: 'SURE-PD3 Investigators, JAMA 2021 (DOI: 10.1001/jama.2021.10207)',
      publishedResult: 'HALTED FOR FUTILITY (No slowing of clinical progression)',
      bottleneck: 'Extracellular urate cannot cross the inner mitochondrial membrane to fix Complex I leakage or stop cytosolic DA-o-Quinone generation under pO₂ = 1.63× hyperoxia.',
      rescueMechanism: 'Instead of scavenging downstream ROS with urate, Stage 0 GBT601 + PT2399 (PDB 5UFP) clamps upstream tissue pO₂ from 1.63× to 0.92× while Stage 1 purges the ROS-producing D mitochondria.',
      gateDay: 27, ambroxol: 0.12, hypoxiaShield: 0.15, kinetin: 0.10, stage2Build: 0.18, creatineMode: 'none',
      trialCurveRegimen: 'levodopa',
      msg: 'SURE-PD3 (Inosine, N=298) failed because raising serum urate acts as a weak downstream antioxidant that cannot prevent midbrain tissue hyperoxia (pO₂ = 1.63×) or purge damaged mtDNA templates.'
    },
    net_pd_ls1: {
      name: '3. NET-PD LS-1 (Creatine Monohydrate 10g/d)',
      nct: 'NCT00449865 · Phase III (N = 1,741, 5 Years)',
      citation: 'Writing Group for the NINDS Exploratory Trials in PD, JAMA 2015 (DOI: 10.1001/jama.2015.120)',
      publishedResult: 'HALTED FOR FUTILITY (N = 1,741; uMtCK oxidized to inactive dimers)',
      bottleneck: 'Discovery 11 (PDB 1QK1): Giving Creatine on Day 1 fails because high Stage-1 mitochondrial ROS/ONOO⁻ oxidizes the catalytic Cys283 residue of octameric uMtCK (leaving only 42.4% active shuttle) and premature ATP buffering suppresses the AMPK/ULK1 energy drop required to trigger mitophagy!',
      rescueMechanism: 'Withholding Creatine during Stage 1 allows AMPK-driven mitophagy to clear D to <8%. Introducing Creatine AFTER the Day-27 Gate preserves 97.0% octameric uMtCK (1QK1) and drives 99.1% TH rescue.',
      gateDay: 14, ambroxol: 0.20, hypoxiaShield: 0.10, kinetin: 0.20, stage2Build: 0.85, creatineMode: 'day1_unphased',
      trialCurveRegimen: 'unstaged_day1',
      msg: 'NET-PD LS-1 (Creatine 10g/d, N=1,741) is the classic example of why Day-1 energy buffering fails: in an unpurged neuron (h=58%), mitochondrial peroxynitrite oxidizes uMtCK Cys283 (PDB 1QK1, dropping octamer activity to 42.4%) while blunting AMPK mitophagy signaling.'
    },
    qe3_coq10: {
      name: '4. QE3 Phase III (High-Dose Coenzyme Q10)',
      nct: 'NCT00740714 · Phase III (N = 600, 16 Months)',
      citation: 'Parkinson Study Group QE3 Investigators, JAMA Neurol 2014',
      publishedResult: 'HALTED FOR FUTILITY (No benefit at 1,200 or 2,400 mg/d)',
      bottleneck: 'Pumping electrons into Complex III via CoQ10 while 58% of mtDNA templates carry Complex I/IV defects and lysosomes are jammed fails to remove the defective D genomes.',
      rescueMechanism: 'CoQ10/Ubiquinol must be deployed in Stage 2 AFTER Stage 1 Kinetin/Urolithin A has purged mutant mtDNA below 8%.',
      gateDay: 14, ambroxol: 0.15, hypoxiaShield: 0.10, kinetin: 0.15, stage2Build: 0.75, creatineMode: 'none',
      trialCurveRegimen: 'unstaged_day1',
      msg: 'QE3 (CoQ10, N=600) failed because supplying electron-transport cofactors to a 58% mutant mitochondrial pool without prior mitophagy clearance leaves the defective genomes intact.'
    },
    fs_zone_pio: {
      name: '5. FS-ZONE (Pioglitazone PPARγ/PGC-1α Agonist)',
      nct: 'NCT01280123 · Phase II (N = 210, 44 Weeks)',
      citation: 'NINDS FS-ZONE Investigators, Lancet Neurol 2015',
      publishedResult: 'FAILED (Unphased Day-1 biogenesis cloned mutant mtDNA!)',
      bottleneck: 'Stimulating PGC-1α mitochondrial biogenesis on Day 1—when 58% of existing mitochondrial templates are mutant (D) with a 1.18× replicative advantage—preferentially clones mutant mitochondria!',
      rescueMechanism: 'Hold all PGC-1α biogenesis agonists until AFTER the Day-27 Biomarker Gate (once h < 8%). Then Stage 2 biogenesis clones 92%+ pure wild-type W mitochondria.',
      gateDay: 14, ambroxol: 0.15, hypoxiaShield: 0.10, kinetin: 0.12, stage2Build: 0.90, creatineMode: 'day1_unphased',
      trialCurveRegimen: 'unstaged_day1',
      msg: 'FS-ZONE (Pioglitazone, N=210) demonstrates the "Unstaged Biogenesis Trap": triggering PGC-1α mitochondrial cloning on Day 1 amplifies the 58% damaged D templates instead of clearing them.'
    },
    adagio_rasagiline: {
      name: '6. ADAGIO (Rasagiline Delayed-Start MAO-B)',
      nct: 'NCT00256204 · Phase III (N = 1,176, 72 Weeks)',
      citation: 'Olanow et al., N Engl J Med 2009',
      publishedResult: '2 mg/d Arm FAILED Delayed-Start Disease Modification',
      bottleneck: 'Inhibiting MAO-B slows enzymatic dopamine breakdown but does not stop non-enzymatic autoxidation of dopamine by excess molecular O₂ (midbrain hyperoxia pO₂ = 1.63×).',
      rescueMechanism: 'Combine low-dose dopamine sparing with Stage 0 GBT601 + PT2399 (PDB 5UFP) tissue O₂ clamping to eliminate non-enzymatic Dopamine-o-Quinone production.',
      gateDay: 27, ambroxol: 0.15, hypoxiaShield: 0.20, kinetin: 0.10, stage2Build: 0.20, creatineMode: 'none',
      trialCurveRegimen: 'levodopa',
      msg: 'ADAGIO (Rasagiline, N=1,176) failed to halt progression because non-enzymatic reaction of cytosolic dopamine with excess midbrain O₂ (pO₂ = 1.63×) continues generating toxic DA-o-Quinone even when MAO-B is blocked.'
    },
    spark_cinpanemab: {
      name: '7. SPARK & PADOVA (Anti-α-Synuclein Antibodies)',
      nct: 'NCT03318523 · Phase II (N = 357, Cinpanemab / Prasinezumab)',
      citation: 'Lang et al., N Engl J Med 2022 (DOI: 10.1056/NEJMoa2209254)',
      publishedResult: 'FAILED Primary Endpoints (Extracellular mAbs miss intracellular GBA1 jam)',
      bottleneck: 'Monoclonal antibodies only bind extracellular α-synuclein in interstitial fluid; they cannot enter the neuron to unclog intracellular GBA1 lysosomes (PDB 2NSX) or stop DA-o-Quinone from cross-linking intracellular α-synuclein.',
      rescueMechanism: 'Intracellular small-molecule chaperone Ambroxol (PDB 2NSX) in Stage 0 restores lysosomal GBA1 activity and clears intracellular α-synuclein oligomers from 1.00× down to 0.19×.',
      gateDay: 27, ambroxol: 0.18, hypoxiaShield: 0.10, kinetin: 0.10, stage2Build: 0.15, creatineMode: 'none',
      trialCurveRegimen: 'levodopa',
      msg: 'SPARK (Cinpanemab, NEJM 2022) failed because extracellular antibodies cannot fix the intracellular root cause of α-synuclein aggregation: clogged GBA1 lysosomes (PDB 2NSX) and intracellular DA-o-Quinone cross-linking.'
    },
    exenatide_pd3: {
      name: '8. Exenatide-PD3 (GLP-1 Receptor Agonist)',
      nct: 'NCT04232969 · Phase III (N = 194, 96 Weeks)',
      citation: 'Vijiaratnam et al., Lancet 2025',
      publishedResult: 'NEGATIVE Phase III (No difference vs. placebo at 96 weeks)',
      bottleneck: 'GLP-1R stimulation modestly dampens microglial inflammation and boosts insulin/PI3K signaling, which actually suppresses ULK1/AMPK autophagy during the window when damaged mitochondria need to be purged.',
      rescueMechanism: 'Reserve anabolic/trophic signaling for Stage 2 AFTER the pulsed fasting + Kinetin Stage 1 purge has cleared mutant mtDNA below 8%.',
      gateDay: 14, ambroxol: 0.25, hypoxiaShield: 0.15, kinetin: 0.18, stage2Build: 0.70, creatineMode: 'none',
      trialCurveRegimen: 'unstaged_day1',
      msg: 'Exenatide-PD3 (Lancet 2025) failed because continuous trophic GLP-1 signaling without a prior Stage 1 mitophagy purge leaves the 58% mutant mtDNA burden intact.'
    },
    nadpark_nr: {
      name: '9. NR-PD / NOPARK (Day-1 Nicotinamide Riboside NAD⁺)',
      nct: 'NCT03816020 · Phase II (Unstaged High-Dose NR 1,000 mg/d)',
      citation: 'Brakedal et al., Cell Metab 2022 / High-Dose NR Follow-up',
      publishedResult: 'Transient metabolic boost, but plateaus without clearing mutant mtDNA',
      bottleneck: 'Giving high-dose NAD⁺ precursors on Day 1 activates SIRT1/PGC-1α biogenesis while h₀ = 58%, feeding both mutant D and wild-type W mitochondria indiscriminately.',
      rescueMechanism: 'Delay high-dose NAD⁺ + SIRT1 (PDB 4ZZJ) + 5-MTHF (PDB 6FCX) until Stage 2 (Day 28+), where it selectively expands the purified (>92% W) mitochondrial pool and drives TET de-aging.',
      gateDay: 14, ambroxol: 0.30, hypoxiaShield: 0.20, kinetin: 0.22, stage2Build: 0.92, creatineMode: 'day1_unphased',
      trialCurveRegimen: 'unstaged_day1',
      msg: 'Day-1 High-Dose NR (NAD⁺ precursor) plateaus clinically because boosting SIRT1/PGC-1α biogenesis before purging mutant mtDNA clones both damaged and healthy mitochondria.'
    },
    biib094_lrrk2: {
      name: '10. BIIB094 / DNL151 (LRRK2 Kinase Inhibitor Monotherapy)',
      nct: 'NCT03976349 · Phase I/II (LRRK2 Kinase Inhibition Alone)',
      citation: 'Jennings et al., Sci Transl Med 2022',
      publishedResult: 'Restores lysosomal Rab phosphorylation, but does not purge existing mutant mtDNA',
      bottleneck: 'Inhibiting hyperactive LRRK2 (PDB 7CK6) takes the foot off the lysosomal brake (Stage 0), but without Stage 1 PINK1 activation (Kinetin) and Stage 2 biogenesis, the existing 58% damaged mitochondria remain in place.',
      rescueMechanism: 'Combine Stage 0 lysosomal unblocking (Ambroxol + LRRK2 normalization) with Stage 1 Kinetin mitophagy and Stage 2 clonal W expansion.',
      gateDay: 27, ambroxol: 0.88, hypoxiaShield: 0.20, kinetin: 0.15, stage2Build: 0.20, creatineMode: 'none',
      trialCurveRegimen: 'no_ambroxol',
      msg: 'LRRK2 monotherapy (PDB 7CK6) solves only Stage 0 (lysosomal trafficking) without executing Stage 1 (selective PINK1 mitophagy) or Stage 2 (clonal wild-type biogenesis).'
    }
  };

  function clip(v, lo, hi) {
    return Math.max(lo, Math.min(hi, v));
  }

  function loadState() {
    try {
      const raw = sessionStorage.getItem(STORAGE_KEY);
      if (raw) {
        const parsed = JSON.parse(raw);
        return Object.assign({}, DEFAULT_STATE, parsed);
      }
    } catch (e) { /* ignore */ }
    return Object.assign({}, DEFAULT_STATE);
  }

  function saveState(st) {
    try {
      const copy = Object.assign({}, st);
      delete copy.trajectories;
      sessionStorage.setItem(STORAGE_KEY, JSON.stringify(copy));
    } catch (e) { /* ignore */ }
  }

  const pdState = loadState();
  pdState.trajectories = {};

  /**
   * 10-State RK4/Euler Biophysical ODE Solver across 90 Days (dt = 0.2d, 451 steps)
   * Supports 5 distinct regimens to answer every mechanistic question:
   *   1. 'levodopa'          — ❌ Standard Symptomatic Care (No Purge, High DA-Quinone)
   *   2. 'unbraked_purge'    — ☠️ Aggressive Purge WITHOUT Stage 0 Shield or Day-27 Brake (Starves & Jams the Neuron!)
   *   3. 'unstaged_day1'     — ⚠️ Unstaged Day-1 Biogenesis + Creatine Stack (Clones Mutant D & Oxidizes uMtCK!)
   *   4. 'no_ambroxol'       — ⚠️ Mitophagy Without Stage 0 GBA1 Lysosomal Unblocking (Autophagosome Traffic Jam!)
   *   5. 'hypoxia_only'      — 💊 GBT601 + PT2399 Oxygen Clamp Alone (Protects DA, But Leaves 52% Mutant mtDNA)
   *   6. 'pd_bcr'            — ✅ Phased PD-BCR Reset (Stage 0 Shield -> Stage 1 Selective Purge + Ketone Bridge -> Day-27 Gate -> Stage 2 Rebuild)
   */
  function solveSingleArmODE(armName) {
    const N = 450;
    const dt = 0.2;

    let W = 1.0 - pdState.h0;
    let D = pdState.h0;
    let H2S = 0.42;
    let Jam = 0.82 * pdState.genoJamScale;
    let pO2 = 1.63 * pdState.genoHyperoxiaScale;
    let DAQ = 0.72;
    let SNCA = 1.00;
    let AKG = 0.62;
    let eAge = 78.0;
    let uMtCK = 48.0;

    const arrDays = new Float32Array(N + 1);
    const arrH = new Float32Array(N + 1);
    const arrW = new Float32Array(N + 1);
    const arrD = new Float32Array(N + 1);
    const arrATP = new Float32Array(N + 1);
    const arrJam = new Float32Array(N + 1);
    const arrSNCA = new Float32Array(N + 1);
    const arrPo2 = new Float32Array(N + 1);
    const arrDAQ = new Float32Array(N + 1);
    const arrAKG = new Float32Array(N + 1);
    const arrEAge = new Float32Array(N + 1);
    const arrUMtCK = new Float32Array(N + 1);
    const arrTH = new Float32Array(N + 1);

    for (let i = 0; i <= N; i++) {
      const day = i * dt;
      arrDays[i] = day;

      const stage = day < pdState.gateDay ? 1 : 2;
      let u_amb = 0.0, u_hyp = 0.0, u_kin = 0.0, u_nad = 0.0, u_epi = 0.0, u_vdac = 0.0;
      let rep_adv = 1.18;
      let ldopa_push = 0.20;
      let cr_stage1_penalty = 0.0;
      let cr_stage2_boost = 0.0;
      let ketone_bridge = 0.0;

      if (armName === 'levodopa') {
        u_amb = 0.05; u_hyp = 0.0; u_kin = 0.04; u_nad = 0.05; u_epi = 0.0; u_vdac = 0.0;
        rep_adv = 1.18; ldopa_push = 0.36;
      } else if (armName === 'unbraked_purge') {
        // Aggressive purge from Day 1 to Day 90 WITHOUT Stage 0 Ambroxol/Hypoxia shield, WITHOUT Ketone bridge, and WITHOUT Day-27 Brake!
        u_amb = 0.08;
        u_hyp = 0.05;
        u_kin = 0.96; // Unbraked continuous mitophagy!
        u_nad = 0.04; // Never switches to Stage 2 biogenesis!
        u_epi = 0.05;
        u_vdac = 0.10;
        ketone_bridge = 0.0;
        rep_adv = 1.15;
        ldopa_push = 0.22;
      } else if (armName === 'unstaged_day1') {
        // Taking NAD+, PGC-1a agonists, and Creatine on Day 1 without purging D first
        u_amb = 0.20;
        u_hyp = 0.12;
        u_kin = 0.15;
        u_nad = 0.90; // Premature Day-1 biogenesis!
        u_epi = 0.25;
        u_vdac = 0.25;
        cr_stage1_penalty = 0.35; // Blunts AMPK mitophagy and oxidizes Cys283
        rep_adv = 1.24; // Mutant D out-replicates W during Day-1 biogenesis!
        ldopa_push = 0.20;
      } else if (armName === 'no_ambroxol') {
        u_amb = 0.06;
        u_hyp = 0.35 * pdState.hypoxiaShield;
        u_kin = stage === 1 ? pdState.kinetin : 0.25 * pdState.kinetin;
        u_nad = stage === 1 ? 0.15 : 0.65 * pdState.stage2Build;
        u_epi = stage === 2 && pdState.deagingEnabled ? 0.45 * pdState.stage2Build : 0.0;
        u_vdac = 0.35;
        ketone_bridge = stage === 1 ? 0.10 : 0.04;
        rep_adv = 1.18; ldopa_push = 0.22;
      } else if (armName === 'hypoxia_only') {
        u_amb = 0.18;
        u_hyp = pdState.hypoxiaShield;
        u_kin = 0.12;
        u_nad = 0.20;
        u_epi = 0.15;
        u_vdac = 0.40;
        rep_adv = 1.14; ldopa_push = 0.15;
      } else {
        // Phased PD-BCR Protocol
        u_amb = pdState.ambroxol;
        u_vdac = pdState.vdacShield;
        if (pdState.creatineMode === 'day1_unphased') {
          cr_stage1_penalty = stage === 1 ? 0.32 : 0.0;
          cr_stage2_boost = stage === 2 ? 0.35 : 0.0;
        } else if (pdState.creatineMode === 'stage2_gated') {
          cr_stage1_penalty = 0.0;
          cr_stage2_boost = stage === 2 ? 1.0 : 0.0;
        }
        if (stage === 1) {
          // Stage 0 (Days 1-14) ramps up Ambroxol & O2 shield; Stage 1 (Days 15-Gate) runs selective mitophagy + Ketone ATP bridge
          const purgeRamp = day < 8.0 ? 0.35 : 1.0;
          u_kin = pdState.kinetin * purgeRamp * (1.0 - 0.55 * cr_stage1_penalty);
          u_nad = 0.12 + 0.45 * cr_stage1_penalty;
          u_hyp = pdState.hypoxiaShield;
          u_epi = pdState.deagingEnabled ? 0.55 * pdState.stage2Build : 0.15;
          ketone_bridge = 0.24 * (0.5 + 0.5 * pdState.ambroxol); // C8-MCT Ketone bridge sustains ATP during Stage 1!
          rep_adv = 1.15;
        } else {
          // Stage 2 (After Gate): Brake mitophagy (-70%), turn ON clonal W biogenesis + Creatine + TET de-aging
          u_kin = 0.30 * pdState.kinetin;
          u_nad = pdState.stage2Build;
          u_hyp = 0.45 * pdState.hypoxiaShield;
          u_epi = pdState.deagingEnabled ? pdState.stage2Build : 0.20;
          ketone_bridge = 0.12;
          rep_adv = 1.03;
        }
        ldopa_push = 0.08;
      }

      // 1. Gut H2S
      const dH2S = (armName === 'pd_bcr' ? 0.035 - 0.38 * H2S : 0.15 - 0.28 * H2S);
      H2S = clip(H2S + dt * dH2S, 0.05, 1.5);

      // 2. Autophagosome-Lysosome Traffic Jam (PDB: 2NSX GBA1)
      const dJam = (0.36 * pdState.genoJamScale * u_kin * (1.0 - 0.80 * u_amb) + 0.15 * pdState.genoJamScale * SNCA * (1.0 - 0.68 * u_amb)) - (0.28 + 0.74 * u_amb) * Jam;
      Jam = clip(Jam + dt * dJam, 0.07, 1.85);
      const lysoFlux = armName === 'unbraked_purge'
        ? 0.52 // Unbraked autophagy chews up mitochondria recklessly even while jamming lysosomes
        : Math.max(0.08, u_amb / (1.0 + 1.35 * Math.max(0.0, Jam - 0.30)));

      // 3. OXPHOS & Tissue Hyperoxia pO2 Mismatch (Mootha & Ichinose Nature Neurosci / PNAS, PDB: 5UFP)
      const tomm40Open = armName === 'pd_bcr' ? (0.35 + 0.60 * u_amb) : 0.25;
      const effOxphos = Math.max(0.12, (W + 0.12 * D) * (0.58 + 0.42 * tomm40Open) * Math.max(0.32, 1.0 - 0.22 * DAQ));
      const targetPo2 = Math.max(0.90, (1.0 + 0.95 * pdState.genoHyperoxiaScale * Math.max(0.0, 0.88 - effOxphos)) * (1.0 - 0.44 * u_hyp));
      pO2 = clip(pO2 + 0.35 * dt * (targetPo2 - pO2), 0.88, 2.05);

      // 4. Cytosolic Dopamine-o-Quinone (DAQ)
      const dDAQ = (0.24 * ldopa_push + 0.32 * Math.max(0.0, pO2 - 1.0) * (1.0 + 0.45 * D)) - (0.30 + 0.52 * u_hyp + 0.25 * u_amb) * DAQ;
      DAQ = clip(DAQ + dt * dDAQ, 0.04, 1.95);

      // 5. Oligomeric alpha-Synuclein (SNCA)
      const dSNCA = (0.16 * Jam + 0.18 * DAQ) * (1.0 - 0.55 * u_amb) - (0.14 + 0.48 * lysoFlux) * SNCA;
      SNCA = clip(SNCA + dt * dSNCA, 0.16, 2.20);

      // 6 & 7. Healthy Wild-Type (W) and Damaged Mutant (D) mtDNA Pools
      const totalMass = W + D;
      const capFactor = Math.max(0.0, 1.0 - totalMass / 1.18);
      const bioDriveW = 0.022 + 0.072 * u_nad * (0.45 + 0.55 * tomm40Open) * (1.0 + 0.14 * cr_stage2_boost);
      const bioDriveD = bioDriveW * rep_adv * (armName === 'unstaged_day1' ? 1.35 : 1.0);

      const purgeD = (0.012 + 0.092 * pdState.genoMitophagyScale * u_kin * lysoFlux) * (1.0 - 0.42 * cr_stage1_penalty);
      // In unbraked purge without Ketone bridge or Day-27 brake, excessive non-selective autophagy also erodes W!
      const purgeW = armName === 'unbraked_purge'
        ? (0.018 + 0.028 * Math.max(0, day - 22) / 68.0)
        : 0.006 * (1.0 + 0.25 * Math.max(0.0, Jam - 0.65));

      W = clip(W + dt * (bioDriveW * W * capFactor - purgeW * W), 0.08, 1.15);
      D = clip(D + dt * (bioDriveD * D * capFactor - purgeD * D), 0.008, 1.10);

      const hPct = (D / (W + D)) * 100.0;

      // 8. Alpha-Ketoglutarate : Succinate Ratio (Krebs Flux for TET1-3)
      const targetAKG = clip(0.42 + 0.78 * W - 0.35 * D + 0.32 * u_epi - 0.18 * H2S, 0.30, 1.65);
      AKG = clip(AKG + 0.28 * dt * (targetAKG - AKG), 0.30, 1.65);

      // 9. Post-Mitotic Neuronal Epigenetic Age (eAge%, Berg et al. 2026 Aging Cell, DOI: 10.1111/acel.70751)
      const tetFlux = u_epi * Math.max(0.0, AKG - 0.72) * (0.55 + 0.45 * u_vdac);
      const eAgeDrift = 0.08 * (D + 0.4 * DAQ);
      eAge = clip(eAge + dt * (eAgeDrift - 0.155 * tetFlux * eAge), 2.4, 98.0);

      // 10. Octameric uMtCK Creatine Shuttle Integrity (PDB: 1QK1 Cys283 Oxidation State)
      const rosLoad = 0.55 * D + 0.45 * Math.max(0.0, pO2 - 1.0) + 0.35 * DAQ;
      const targetUMtCK = clip(
        98.5 - 52.0 * rosLoad - 14.0 * cr_stage1_penalty + 14.5 * cr_stage2_boost * W,
        22.0,
        98.5
      );
      uMtCK = clip(uMtCK + 0.30 * dt * (targetUMtCK - uMtCK), 22.0, 98.5);

      // Instantaneous Neuronal Pacemaker ATP (Survival Floor = 0.45x)
      const crShuttleBonus = (uMtCK / 100.0) * (0.16 * cr_stage2_boost + 0.06);
      const instATP = clip(
        0.88 * effOxphos + ketone_bridge + crShuttleBonus - 0.16 * Math.max(0.0, Jam - 0.35) - 0.12 * Math.max(0.0, DAQ - 0.40),
        0.18,
        1.52
      );

      // Composite Dopaminergic TH & Motor Rescue Score (%)
      const survivalPenalty = instATP < 0.45 ? (0.45 - instATP) * 140.0 : 0.0;
      const thScore = clip(
        18.0
        + 54.0 * W
        - 32.0 * D
        - 16.0 * Math.max(0.0, Jam - 0.22)
        - 18.0 * Math.max(0.0, pO2 - 1.02)
        + 12.0 * (1.0 - eAge / 100.0)
        + 9.5 * (uMtCK / 100.0) * cr_stage2_boost
        - survivalPenalty,
        4.0,
        99.4
      );

      arrH[i] = hPct;
      arrW[i] = W;
      arrD[i] = D;
      arrATP[i] = instATP;
      arrJam[i] = Jam;
      arrSNCA[i] = SNCA;
      arrPo2[i] = pO2;
      arrDAQ[i] = DAQ;
      arrAKG[i] = AKG;
      arrEAge[i] = eAge;
      arrUMtCK[i] = uMtCK;
      arrTH[i] = thScore;
    }

    return {
      arrDays, arrH, arrW, arrD, arrATP, arrJam, arrSNCA, arrPo2, arrDAQ, arrAKG, arrEAge, arrUMtCK, arrTH
    };
  }

  function solveAllArms() {
    ['levodopa', 'unbraked_purge', 'unstaged_day1', 'no_ambroxol', 'hypoxia_only', 'pd_bcr'].forEach(arm => {
      pdState.trajectories[arm] = solveSingleArmODE(arm);
    });
    saveState(pdState);
  }

  function applyGenotypeProfile(key) {
    const g = GENOTYPE_PROFILES[key];
    if (!g) return;
    pdState.genotype = key;
    pdState.h0 = g.h0;
    pdState.gateDay = g.gateDay;
    pdState.genoJamScale = g.genoJamScale;
    pdState.genoHyperoxiaScale = g.genoHyperoxiaScale;
    pdState.genoMitophagyScale = g.genoMitophagyScale;
    solveAllArms();
  }

  /**
   * Renders the unified Top Header & 4-Page Parkinson's Academic Suite Bar across all PD pages
   */
  function renderPdSuiteHeader(activePageNum) {
    const container = document.getElementById('pdSuiteHeaderContainer');
    if (!container) return;

    const pages = [
      {
        num: 1,
        href: 'parkinsons_simulator.html',
        badge: 'PAGE 1 · CORE MODEL & SURVIVAL',
        title: '1. Pathophysiology & 4-Phase Protocol',
        sub: 'Selective mitophagy vs. starvation floor · 4D Manifold'
      },
      {
        num: 2,
        href: 'parkinsons_patients.html',
        badge: 'PAGE 2 · GENETICS & DOSING',
        title: '2. Patient Genetics & 10-State ODE',
        sub: '5 Subtypes · Custom Labs · 10 Live State Curves'
      },
      {
        num: 3,
        href: 'parkinsons_trials.html',
        badge: 'PAGE 3 · 10-TRIAL AUDIT',
        title: '3. Audit of 10 Failed PD Trials',
        sub: 'Why STEADY-PD III, SURE-PD3 & NET-PD Creatine failed'
      },
      {
        num: 4,
        href: 'parkinsons_breakthroughs.html',
        badge: 'PAGE 4 · 2025–2026 & 3D PDBs',
        title: '4. Breakthroughs, Sweeps & 3D Atlas',
        sub: 'Hyperoxia (5UFP) · Berg De-Aging · uMtCK (1QK1)'
      },
      {
        num: 5,
        href: 'parkinsons_cell_ai_translation.html',
        badge: 'PAGE 5 · CELL (2026) TRANSLATION',
        title: '5. Blood Biosensor, 4-Cell Niche & Trial',
        sub: 'Liquid Biopsy · Astrocyte Donation · 5-Organ Safety · N=500 Trial'
      }
    ];

    container.innerHTML = `
      <!-- Global Platform Navigation Bar -->
      <header class="border-b border-slate-300 bg-[#f6f3ec]/95 backdrop-blur sticky top-0 z-40 shadow-sm">
        <div class="max-w-7xl mx-auto px-4 sm:px-6 py-2.5 flex flex-wrap items-center justify-between gap-3">
          <div class="flex items-center gap-3">
            <a href="../index.html" class="flex items-center gap-2 font-extrabold text-base tracking-tight text-slate-900 hover:text-emerald-800 transition">
              <span class="w-2.5 h-2.5 rounded-full bg-emerald-600"></span>
              <span>MitoRestorer v3.0</span>
            </a>
            <span class="text-slate-300">|</span>
            <span class="text-xs font-mono font-bold px-2.5 py-0.5 rounded bg-emerald-950 text-emerald-200 border border-emerald-500/40">
              🧠 FLAGSHIP #1: PARKINSON'S ACADEMIC SUITE (PD-BCR · 10-STATE ODE)
            </span>
          </div>

          <!-- Unified Top-Level Platform Switcher (Resolves Cohort vs. Dedicated Page Confusion) -->
          <div class="flex flex-wrap items-center gap-2 text-xs font-semibold">
            <a href="../index.html" class="px-2.5 py-1.5 rounded-lg bg-white border border-slate-300 text-slate-800 hover:border-emerald-600 transition">
              🏠 Platform Overview
            </a>
            <a href="apis_honeybee_simulator.html" class="px-2.5 py-1.5 rounded-lg bg-amber-50 border border-amber-300 text-amber-900 font-bold hover:bg-amber-100 transition">
              🐝 Flagship #2: Honeybee &amp; Coral NPU
            </a>
            <details class="relative group">
              <summary class="list-none cursor-pointer px-2.5 py-1.5 rounded-lg bg-white border border-slate-300 text-slate-800 hover:border-cyan-600 transition font-semibold">
                🧬 Other Clinical &amp; Bioreactor Cohorts (5) ▾
              </summary>
              <div class="absolute right-0 mt-1.5 w-72 rounded-xl bg-white border border-slate-300 shadow-2xl p-2 z-50 space-y-1 text-xs">
                <div class="px-2 py-1 text-[10px] font-mono font-bold text-slate-500 uppercase border-b border-slate-200">
                  Dedicated Single-Cohort Workstations
                </div>
                <a href="cart_bioreactor_simulator.html" class="block px-2.5 py-1.5 rounded-lg hover:bg-slate-100 text-slate-900 font-semibold">
                  🧫 1. Ex Vivo CAR-T Bioreactor (240h GMP)
                </a>
                <a href="alzheimers_simulator.html" class="block px-2.5 py-1.5 rounded-lg hover:bg-slate-100 text-slate-900">
                  🧠 2. Prodromal APOE4 Alzheimer's Disease
                </a>
                <a href="epilepsy_simulator.html" class="block px-2.5 py-1.5 rounded-lg hover:bg-slate-100 text-slate-900">
                  ⚡ 3. Refractory Epilepsy (PV+ Interneurons)
                </a>
                <a href="depression_simulator.html" class="block px-2.5 py-1.5 rounded-lg hover:bg-slate-100 text-slate-900">
                  🔬 4. Treatment-Resistant Depression (MDD)
                </a>
                <a href="sarcopenia_simulator.html" class="block px-2.5 py-1.5 rounded-lg hover:bg-slate-100 text-slate-900">
                  💪 5. Sarcopenia &amp; Metabolic Syndrome
                </a>
                <div class="border-t border-slate-200 pt-1 mt-1">
                  <a href="bcr_simulator.html" class="block px-2.5 py-1.5 rounded-lg hover:bg-slate-100 text-emerald-800 font-bold">
                    📐 Cross-Cohort Comparative Master View
                  </a>
                </div>
              </div>
            </details>
            <a href="science_radar_sandbox.html" class="px-2.5 py-1.5 rounded-lg bg-white border border-slate-300 text-slate-800 hover:border-cyan-600 transition">
              🛰️ AI Science Radar
            </a>
            <a href="manuscript_reader.html" class="px-2.5 py-1.5 rounded-lg bg-emerald-800 text-white font-bold hover:bg-emerald-700 transition">
              📄 Read Academic Preprint (32 DOIs)
            </a>
          </div>
        </div>

        <!-- 5-Page Parkinson's Academic Module Stepper Bar -->
        <div class="border-t border-slate-300/80 bg-[#eee9dd]">
          <div class="max-w-7xl mx-auto px-4 sm:px-6 py-2">
            <div class="grid grid-cols-1 sm:grid-cols-2 lg:grid-cols-5 gap-2">
              ${pages.map(p => {
                const isActive = p.num === activePageNum;
                return `
                  <a href="${p.href}" class="block p-2.5 rounded-xl border transition ${
                    isActive
                      ? 'bg-white border-2 border-emerald-700 shadow-sm'
                      : 'bg-[#f8f6f0] border-slate-300 hover:bg-white hover:border-emerald-600'
                  }">
                    <div class="flex items-center justify-between gap-1">
                      <span class="font-mono text-[10px] font-bold ${isActive ? 'text-emerald-800' : 'text-slate-600'}">${p.badge}</span>
                      ${isActive ? '<span class="px-1.5 py-0.5 rounded bg-emerald-800 text-white font-mono text-[9px] font-bold">ACTIVE PAGE</span>' : '<span class="text-[10px] font-mono text-emerald-700 font-bold">Open →</span>'}
                    </div>
                    <div class="font-extrabold text-xs sm:text-[13px] text-slate-900 mt-0.5">${p.title}</div>
                    <div class="text-[11px] text-slate-600 mt-0.5 leading-snug">${p.sub}</div>
                  </a>
                `;
              }).join('')}
            </div>
          </div>
        </div>
      </header>
    `;
  }

  global.PDBCR = {
    pdState,
    GENOTYPE_PROFILES,
    PHASE_DEFINITIONS,
    TRIAL_AUDIT_PRESETS,
    clip,
    solveSingleArmODE,
    solveAllArms,
    applyGenotypeProfile,
    saveState,
    renderPdSuiteHeader
  };
})(window);

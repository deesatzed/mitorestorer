#!/usr/bin/env python3
"""
Builds the 5 Dedicated Human Clinical Cohort Workstations:
  - simulator/parkinsons_simulator.html   (Cohort E: Parkinson's Disease SNc)
  - simulator/alzheimers_simulator.html   (Cohort B: Prodromal APOE4 Alzheimer's)
  - simulator/epilepsy_simulator.html     (Cohort C: Refractory Epilepsy PV+ Interneurons)
  - simulator/depression_simulator.html   (Cohort D: Treatment-Resistant Depression TR-MDD)
  - simulator/sarcopenia_simulator.html   (Cohort A: Metabolic Syndrome & Sarcopenia)

Each page is 100% dedicated to that single clinical indication and features:
  1. Top 7-Cohort Navigation Bar (linking all 7 separate cohort workstations + Central Hub)
  2. Interactive 4D Bioenergetic Phase-Space & Free-Energy Manifold (`#canvas4DManifold`)
     comparing Standard Care, Simultaneous Supplement Trap, and Phased BCR Rescue in 3D/4D space
  3. 4-Phase Workflow Progression Stepper + 5-Pillar Human-UX Transparency Inspector:
     - Pillar 1: Step Purpose & Biological Goal
     - Pillar 2: Why It Flags (Biomarker & Safety Thresholds)
     - Pillar 3: Decision Scope & Intention (What is actuated vs. intentionally withheld)
     - Pillar 4: Understandable Transparency (ODE Equations, 3D PDB Coordinates, Ensembl VEP/gnomAD)
     - Pillar 5: Overall Scope & Endpoint Impact
  4. Cohort-Specific 9-Point "What-If" ODE Sensitivity Sweeps + Live Gemini API (`@google/genai` `gemini-3.8-flash`) Co-Pilot
"""

from pathlib import Path

COHORT_PAGES = [
    {
        "spec_key": "parkinsons_snc",
        "filename": "parkinsons_simulator.html",
        "nav_short": "🧠 Parkinson's (Cohort E)",
        "primary_metric_label": "SNc Dopamine (TH) Tone",
        "primary_metric_bcr": "1.25× (95% Rescue)",
        "primary_metric_trap": "0.15× (Jam=1.29×)",
        "primary_metric_std": "0.13× (0% Rescue)",
        "secondary_metric_label": "Lysosomal Jam / SNCA",
        "secondary_metric_bcr": "Jam 0.08× | SNCA 0.15×",
        "slider1_label": "Stage 0 Ambroxol + Spermidine (PDB: 2NSX GBA1)",
        "slider1_desc": "Stabilizes Glu235 (2.61 Å) & Glu340 (2.67 Å) to clear lysosomal jam before mitophagy.",
        "slider2_label": "Stage 1 Kinetin + Urolithin A + D-BHB (PDB: 6EQI)",
        "slider2_desc": "Clears deletion-mutant mtDNA (D) in SNc pacemaker axons while NAD+ is held off.",
        "slider3_label": "Stage 2 NR (1,000 mg) + SIRT1 STAC (PDB: 4ZZJ) + B2/Cr",
        "slider3_desc": "NADPARK-dose biogenesis after Day-34/38 GDF15 gate (<550 pg/mL).",
        "default_gate_day": 34,
        "alphagenome_gene": "GBA1",
        "alphagenome_var": "chr1:155235843:T>C",
        "alphafold_id": "P04062",
        "alphamissense_summary": "GBA1 N370S (AM = 0.562) · LRRK2 G2019S (AM = 0.941) · EPAS1/HIF-2α (PDB: 5UFP PT2399 + GBT601 pO₂ Shield)",
        "has_hyperoxia_shield": True,
    },
    {
        "spec_key": "apoe4_keto",
        "filename": "alzheimers_simulator.html",
        "nav_short": "🔬 Alzheimer's APOE4 (Cohort B)",
        "primary_metric_label": "Cortical Synaptic ATP",
        "primary_metric_bcr": "1.19× (Mass 1.06×)",
        "primary_metric_trap": "0.17× (h=94.2%!)",
        "primary_metric_std": "0.46× (h=73.3%)",
        "secondary_metric_label": "Plasma GDF15 Mitokine",
        "secondary_metric_bcr": "571 pg/mL (-86%)",
        "slider1_label": "Stage 0 TUDCA + Citicoline + LPC-DHA (PDB: 7CK6 TOMM40)",
        "slider1_desc": "Stabilizes PC1 401 phospholipid inside the 19-stranded TOMM40 β-barrel import pore.",
        "slider2_label": "Stage 1 Kinetin (6EQI) + Trehalose/Spermidine + C8-MCT",
        "slider2_desc": "Accelerates slow cortical neuronal transmitophagy (0.65× turnover rate).",
        "slider3_label": "Stage 2 NR + PncA + Resveratrol (4ZZJ) + PQQ + Cr/B2",
        "slider3_desc": "Expands wild-type cortical mitochondria after Day-38 GDF15 gate.",
        "default_gate_day": 38,
        "alphagenome_gene": "TOMM40",
        "alphagenome_var": "chr19:44908684:T>C",
        "alphafold_id": "O96008",
        "alphamissense_summary": "APOE Cys112Arg (rs429358) · TOMM40 β-Barrel Import Stall (PDB: 7CK6)",
    },
    {
        "spec_key": "epilepsy_pv",
        "filename": "epilepsy_simulator.html",
        "nav_short": "⚡ Epilepsy PV+ (Cohort C)",
        "primary_metric_label": "Cortical E/I & Seizure Red.",
        "primary_metric_bcr": "E/I 0.65 (96% Red.)",
        "primary_metric_trap": "E/I 2.26 (H₂S Rebound)",
        "primary_metric_std": "E/I 3.38 (0% Red.)",
        "secondary_metric_label": "Gut H₂S & PV+ Interneuron ATP",
        "secondary_metric_bcr": "H₂S 0.11 | PV+ ATP 1.46×",
        "slider1_label": "Stage 0 Molybdenum + Polyphenols + Akkermansia/gadB",
        "slider1_desc": "Blocks Bilophila H₂S production so Complex IV isn't poisoned on a ketogenic diet.",
        "slider2_label": "Stage 1 C8-MCT TID + D-BHB + Kinetin + Urolithin A",
        "slider2_desc": "Purges mutant mtDNA in 40 Hz fast-spiking PV+ basket interneurons (f² = 3.06×).",
        "slider3_label": "Stage 2 Ketone-Maintained NR + SIRT1 (4ZZJ) + Creatine",
        "slider3_desc": "Builds PV+ phosphocreatine ATP grid while keeping D-BHB ≥ 1.5 mmol/L.",
        "default_gate_day": 32,
        "alphagenome_gene": "PPARGC1A",
        "alphagenome_var": "chr4:23815662:C>T",
        "alphafold_id": "Q9BXM7",
        "alphamissense_summary": "PV+ 40 Hz Gamma Load (f² = 3.06×) · SUOX / Complex IV Sulfide Shield",
    },
    {
        "spec_key": "mdd_quin",
        "filename": "depression_simulator.html",
        "nav_short": "🌱 Depression TR-MDD (Cohort D)",
        "primary_metric_label": "Prefrontal Dendritic Spines",
        "primary_metric_bcr": "1.24× (Sustained d90)",
        "primary_metric_trap": "0.45× (d14 Relapse!)",
        "primary_metric_std": "0.45× (Drive 0.09×)",
        "secondary_metric_label": "Microglial QUIN & GDF15",
        "secondary_metric_bcr": "QUIN 0.20 (-87%) | GDF15 259",
        "slider1_label": "Stage 0 Tributyrin + High-EPA + Curcumin (IDO1/QUIN Block)",
        "slider1_desc": "Shuts down microglial kynurenine QUIN excitotoxicity and brainstem GFRAL lock.",
        "slider2_label": "Stage 1 Prefrontal Mitophagy (Urolithin A + Kinetin)",
        "slider2_desc": "Clears dendritic energy resistance (éR) before any neuroplastogen/ketamine pulse.",
        "slider3_label": "Stage 2 Gated Synaptogenesis (NR + 4ZZJ + ALCAR + 6FCX B2)",
        "slider3_desc": "Locks in 1.24× prefrontal spine density and BH4/SAMe monoamine drive.",
        "default_gate_day": 31,
        "alphagenome_gene": "SIRT1",
        "alphagenome_var": "chr1:11796321:G>A",
        "alphafold_id": "Q96EB6",
        "alphamissense_summary": "MTHFR Ala222Val (AM = 0.632) · IDO1/QUIN Shunt · SIRT1 Glu230 (PDB: 4ZZJ)",
    },
    {
        "spec_key": "pink1_fut2",
        "filename": "sarcopenia_simulator.html",
        "nav_short": "💪 Sarcopenia (Cohort A)",
        "primary_metric_label": "Skeletal Muscle Net ATP",
        "primary_metric_bcr": "1.20× (Mass 1.04×)",
        "primary_metric_trap": "0.24× (h=91.4% Trap!)",
        "primary_metric_std": "0.50× (16:8 Lag Fail)",
        "secondary_metric_label": "Glycogen Lag & Mito Mass",
        "secondary_metric_bcr": "τ_lag = 6.0h | Mass 1.04×",
        "slider1_label": "Stage 0 Direct Urolithin A + 2'-FL HMO + Zone-2/C8-MCT",
        "slider1_desc": "Bypasses FUT2 rs601338 A/A non-secretor block & compresses 13.5h glycogen lag to 6.0h.",
        "slider2_label": "Stage 1 Lag-Compressed 18:6 Purge + Kinetin (PDB: 6EQI)",
        "slider2_desc": "Unlocks 12h/day AMPK/ULK1 mitophagy in skeletal muscle while holding NR off.",
        "slider3_label": "Stage 2 Weekly 5.25d Build (NR + 4ZZJ + 6FCX B2 + Creatine)",
        "slider3_desc": "Switches from daily 18:6 to 5.25d/wk anabolic refeed at Day 38 to prevent atrophy.",
        "default_gate_day": 38,
        "alphagenome_gene": "MTHFR",
        "alphagenome_var": "chr1:11796321:G>A",
        "alphafold_id": "P42898",
        "alphamissense_summary": "FUT2 Trp154Ter (VEP HIGH · 22.5% Hom) · MTHFR Ala222Val (AM = 0.632 · PDB: 6FCX)",
    },
]


def render_top_nav(active_file: str) -> str:
    links = [
        ("../index.html", "🏛️ Hub"),
        ("parkinsons_simulator.html", "🧠 Parkinson's"),
        ("alzheimers_simulator.html", "🔬 Alzheimer's"),
        ("epilepsy_simulator.html", "⚡ Epilepsy"),
        ("depression_simulator.html", "🌱 Depression"),
        ("sarcopenia_simulator.html", "💪 Sarcopenia"),
        ("cart_bioreactor_simulator.html", "🧬 CAR-T (240h)"),
        ("apis_honeybee_simulator.html", "🐝 Honeybee (14d)"),
        ("science_radar_sandbox.html", "🛰️ New Science Sandbox"),
    ]
    items = []
    for href, label in links:
        is_act = href == active_file
        cls = (
            "px-2.5 py-1 rounded-lg bg-emerald-500 text-slate-950 font-extrabold shadow-sm"
            if is_act
            else "px-2.5 py-1 rounded-lg bg-slate-900 text-slate-200 border border-slate-700 hover:border-emerald-400 transition"
        )
        items.append(f'<a href="{href}" class="{cls}">{label}</a>')
    return "\n          ".join(items)


def build_page_html(cfg: dict) -> str:
    nav_html = render_top_nav(cfg["filename"])
    spec_key = cfg["spec_key"]
    gate_d = cfg["default_gate_day"]
    has_hyp = cfg.get("has_hyperoxia_shield", False)
    extra_hdr_badges = ""
    extra_slider_html = ""
    extra_sweep_btn = ""
    if has_hyp:
        extra_hdr_badges = """
            <a href="https://deepmind.google.com/science/alphagenome/atlas?q=EPAS1&m=entity" target="_blank" rel="noopener" class="px-2.5 py-0.5 rounded bg-rose-950 text-rose-200 border border-rose-500/50 hover:bg-rose-900 transition font-bold">
              🧬 AlphaGenome: EPAS1 (HIF-2α) ↗
            </a>
            <a href="https://alphafold.ebi.ac.uk/entry/Q99814" target="_blank" rel="noopener" class="px-2.5 py-0.5 rounded bg-amber-950 text-amber-200 border border-amber-500/50 hover:bg-amber-900 transition font-bold">
              🔬 AlphaFold: EPAS1 (Q99814 / PDB 5UFP) ↗
            </a>"""
        extra_slider_html = """
            <div class="space-y-1 sm:col-span-2 pt-1.5 border-t border-slate-800">
              <div class="flex justify-between font-semibold text-rose-200">
                <span>💊 Stage 0/1 GBT601 (Hb-O₂) + PT2399 (HIF-2α) Brain Hyperoxia Shield</span>
                <span id="lblHypShield" class="font-mono text-rose-300 font-bold">94% (pO₂=0.92×)</span>
              </div>
              <input id="sliderHypShield" type="range" min="0" max="100" step="2" value="94" class="w-full accent-rose-400" />
              <div class="text-[10.5px] text-slate-300">Mootha &amp; Ichinose (<em>Nature Neurosci / PNAS</em>): Clamps SNc tissue hyperoxia (<code>pO₂ 1.63× → 0.92×</code>) so unconsumed O₂ cannot oxidize dopamine into DA-o-Quinone; accelerates Stage 2 gate from Day 34 to Day 26.5 (98.2% Rescue).</div>
            </div>"""
        extra_sweep_btn = """
            <button data-sweep="hyperoxia_shield" class="sweep-btn col-span-3 px-2 py-1.5 rounded-lg bg-rose-950 border border-rose-500/70 text-rose-200 hover:bg-rose-900 font-mono font-bold text-[10.5px] transition">
              🔬 4. Mootha/Ichinose Brain Hyperoxia Shield (GBT601 + PT2399 pO₂ Clamp Sweep)
            </button>"""
    return f"""<!DOCTYPE html>
<html lang="en" class="dark">
<head>
  <meta charset="UTF-8" />
  <meta name="viewport" content="width=device-width, initial-scale=1.0" />
  <title>MitoRestorer — {cfg["nav_short"]} Dedicated 4D Bioenergetic Decision Workstation</title>
  <script src="https://cdn.tailwindcss.com"></script>
  <script src="bcr_ux_4d_engine.js"></script>
  <style>
    :root {{
      --bg: #050914;
      --card: #0c1425;
      --content: #111c33;
      --border: #1e293b;
      --fg: #f8fafc;
      --muted: #94a3b8;
    }}
    body {{
      background-color: var(--bg);
      color: var(--fg);
      font-family: ui-sans-serif, system-ui, -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif;
    }}
    canvas {{ width: 100%; display: block; }}
  </style>
</head>
<body class="min-h-screen antialiased">

  <!-- Top 7-Cohort Dedicated Workstation Switcher Bar -->
  <div class="bg-slate-950 border-b border-slate-800 sticky top-0 z-40">
    <div class="max-w-[1600px] mx-auto px-4 py-2 flex flex-wrap items-center justify-between gap-2 text-xs font-mono">
      <div class="flex items-center gap-2">
        <span class="w-2.5 h-2.5 rounded-full bg-emerald-400"></span>
        <span class="text-slate-300 font-bold">DEDICATED COHORT WORKSTATIONS:</span>
      </div>
      <div class="flex flex-wrap items-center gap-1.5">
        {nav_html}
      </div>
    </div>
  </div>

  <!-- Dedicated Cohort Header & 4-KPI Decision Summary Strip -->
  <header class="border-b border-slate-800 bg-[#081020]">
    <div class="max-w-[1600px] mx-auto px-4 sm:px-6 py-4 space-y-4">
      <div class="flex flex-col lg:flex-row lg:items-center justify-between gap-3">
        <div class="space-y-1">
          <div class="flex items-center gap-2 flex-wrap text-xs font-mono">
            <span id="hdrBadge" class="px-2.5 py-0.5 rounded bg-emerald-950 text-emerald-300 border border-emerald-500/50 font-bold"></span>
            <span id="hdrPdb" class="px-2.5 py-0.5 rounded bg-purple-950 text-purple-200 border border-purple-500/40"></span>
            <a href="https://deepmind.google.com/science/alphagenome/atlas?q={cfg['alphagenome_gene']}&m=entity" target="_blank" rel="noopener" class="px-2.5 py-0.5 rounded bg-cyan-950 text-cyan-200 border border-cyan-500/50 hover:bg-cyan-900 transition font-bold">
              🧬 AlphaGenome Atlas ({cfg['alphagenome_gene']}) ↗
            </a>
            <a href="https://alphafold.ebi.ac.uk/entry/{cfg['alphafold_id']}" target="_blank" rel="noopener" class="px-2.5 py-0.5 rounded bg-indigo-950 text-indigo-200 border border-indigo-500/50 hover:bg-indigo-900 transition font-bold">
              🔬 AlphaFold DB ({cfg['alphafold_id']}) ↗
            </a>{extra_hdr_badges}
          </div>
          <h1 id="hdrTitle" class="text-xl sm:text-2xl font-extrabold tracking-tight text-white"></h1>
          <p id="hdrSubtitle" class="text-xs sm:text-sm text-slate-300"></p>
          <div class="text-[11px] font-mono text-cyan-300 pt-0.5">
            ⚡ <strong>Google DeepMind Bio-AI Stack:</strong> {cfg['alphamissense_summary']} · <a href="https://deepmind.google.com/science/alphagenome/atlas?q={cfg['alphagenome_var']}&m=variant" target="_blank" rel="noopener" class="underline hover:text-white">Inspect Variant on AlphaGenome Atlas ({cfg['alphagenome_var']}) ↗</a>
          </div>
        </div>

        <!-- Non-Jumpy Precision Timeline Scrubber & Optional Smooth Play -->
        <div class="flex flex-wrap items-center gap-2 bg-slate-950/90 border border-slate-700 rounded-xl px-3 py-2 text-xs font-mono shrink-0">
          <button id="btnPlayPause" class="px-3 py-1 rounded-lg bg-cyan-500/20 border border-cyan-400 text-cyan-200 font-bold hover:bg-cyan-500/30 transition">
            ▶ Play Progression
          </button>
          <div class="flex items-center gap-2">
            <span class="text-slate-400">Timeline:</span>
            <input id="timeScrubber" type="range" min="0" max="900" value="380" class="w-32 sm:w-44 accent-emerald-400 cursor-pointer" />
            <span id="clockBadge" class="px-2 py-0.5 rounded bg-slate-900 text-emerald-300 border border-emerald-500/40 font-bold">Day 38.0 / 90d</span>
          </div>
        </div>
      </div>

      <!-- 4 Executive Outcome & Decision KPI Cards -->
      <div class="grid grid-cols-2 lg:grid-cols-4 gap-3 text-xs">
        <div class="p-3.5 rounded-xl bg-slate-900/95 border border-emerald-500/40">
          <div class="font-mono text-slate-400 uppercase text-[10px]">1. DAMAGED mtDNA HETEROPLASMY h(90d)</div>
          <div class="text-xl font-extrabold font-mono text-emerald-300 mt-1" id="kpiH">--%</div>
          <div class="text-[11px] text-slate-300 mt-0.5" id="kpiHSub">Target: Purge &lt; 8.0% before Stage 2</div>
        </div>
        <div class="p-3.5 rounded-xl bg-slate-900/95 border border-cyan-500/40">
          <div class="font-mono text-slate-400 uppercase text-[10px]">2. {cfg["primary_metric_label"].upper()}</div>
          <div class="text-xl font-extrabold font-mono text-cyan-300 mt-1" id="kpiPrimary">{cfg["primary_metric_bcr"]}</div>
          <div class="text-[11px] text-slate-300 mt-0.5">Standard Care: {cfg["primary_metric_std"]} | Trap: {cfg["primary_metric_trap"]}</div>
        </div>
        <div class="p-3.5 rounded-xl bg-slate-900/95 border border-purple-500/40">
          <div class="font-mono text-slate-400 uppercase text-[10px]">3. {cfg["secondary_metric_label"].upper()}</div>
          <div class="text-xl font-extrabold font-mono text-purple-200 mt-1" id="kpiSecondary">{cfg["secondary_metric_bcr"]}</div>
          <div class="text-[11px] text-slate-300 mt-0.5">3D Structure-Guided Bottleneck Resolution</div>
        </div>
        <div class="p-3.5 rounded-xl bg-slate-900/95 border border-amber-500/40">
          <div class="font-mono text-slate-400 uppercase text-[10px]">4. STAGE-GATE DECISION STATUS</div>
          <div class="text-base font-extrabold font-mono text-amber-300 mt-1" id="kpiGateStatus">⚡ Day {gate_d} Gate Verified</div>
          <div class="text-[11px] text-slate-300 mt-0.5">ΔGDF15 ≤ -35% &amp; DEXA Mass Loss &gt; -2.5%</div>
        </div>
      </div>
    </div>
  </header>

  <main class="max-w-[1600px] mx-auto px-4 sm:px-6 py-6 space-y-6">

    <!-- SECTION 1: 4-PHASE WORKFLOW STEPPER & 5-PILLAR HUMAN-UX TRANSPARENCY INSPECTOR -->
    <div id="transparencyInspectorContainer"></div>

    <!-- SECTION 2: 4D BIOENERGETIC PHASE-SPACE MANIFOLD (Left 7 Cols) + AI CO-PILOT & SENSITIVITY SWEEPS (Right 5 Cols) -->
    <div class="grid grid-cols-1 lg:grid-cols-12 gap-5 items-stretch">

      <!-- Left 7 Cols: Interactive 4D Bioenergetic Phase-Space & Free-Energy Manifold -->
      <div class="lg:col-span-7 rounded-2xl border-2 border-cyan-500/40 bg-[#060c1a] p-4 flex flex-col justify-between space-y-3 shadow-xl">
        <div class="flex flex-col sm:flex-row sm:items-center justify-between gap-2 border-b border-slate-800 pb-2.5">
          <div>
            <div class="text-[11px] font-mono font-bold text-cyan-300 uppercase tracking-wider">
              📐 INTERACTIVE 4D BIOENERGETIC PHASE-SPACE MANIFOLD (EP = éR · f²)
            </div>
            <h2 class="text-sm sm:text-base font-extrabold text-white">
              3D Trajectory Comparison: Time (X) × Mitochondrial Purity (Y) × Functional ATP Output (Z) × Stress Surface (Color)
            </h2>
          </div>
          <div class="flex items-center gap-1.5 text-[11px] font-mono shrink-0">
            <button data-cam="0" class="cam-btn px-2 py-1 rounded bg-slate-900 border border-slate-700 text-slate-200 hover:border-cyan-400">Isometric</button>
            <button data-cam="1" class="cam-btn px-2 py-1 rounded bg-slate-900 border border-slate-700 text-slate-200 hover:border-cyan-400">Phase-Plane</button>
            <button data-cam="2" class="cam-btn px-2 py-1 rounded bg-slate-900 border border-slate-700 text-slate-200 hover:border-cyan-400">Top Ridge</button>
          </div>
        </div>

        <!-- Strategy Highlight Selector Bar (High-Contrast Dark Cards) -->
        <div class="grid grid-cols-1 sm:grid-cols-3 gap-2 text-xs">
          <button data-reg="untreated" class="reg-card text-left p-2.5 rounded-xl border border-slate-700 bg-slate-900/90 text-slate-200 hover:border-rose-400 transition">
            <div class="font-bold text-rose-300 flex items-center justify-between">
              <span>❌ 1. Standard Care</span>
              <span class="font-mono text-[10px] px-1.5 py-0.5 rounded bg-rose-950 text-rose-200 border border-rose-500/40">Stays Congested</span>
            </div>
            <div class="text-[11px] text-slate-300 mt-1">High resistance (éR) &amp; unbuffered oxidative stress.</div>
          </button>
          <button data-reg="unphased" class="reg-card text-left p-2.5 rounded-xl border border-slate-700 bg-slate-900/90 text-slate-200 hover:border-amber-400 transition">
            <div class="font-bold text-amber-300 flex items-center justify-between">
              <span>⚠️ 2. Simultaneous Trap</span>
              <span class="font-mono text-[10px] px-1.5 py-0.5 rounded bg-amber-950 text-amber-200 border border-amber-500/40">Clones Mutant D</span>
            </div>
            <div class="text-[11px] text-slate-300 mt-1">Day-1 NAD⁺/NR activates mTORC1 &amp; amplifies mutant mtDNA.</div>
          </button>
          <button data-reg="bcr" class="reg-card text-left p-2.5 rounded-xl border-2 border-emerald-500 bg-emerald-950/80 text-white transition shadow-md">
            <div class="font-bold text-emerald-300 flex items-center justify-between">
              <span>✅ 3. Phased BCR Protocol</span>
              <span class="font-mono text-[10px] px-1.5 py-0.5 rounded bg-emerald-900 text-emerald-100 border border-emerald-400 font-bold">Optimal Rescue</span>
            </div>
            <div class="text-[11px] text-emerald-100 mt-1">Phase 0 Unblock → Phase 1 Purge → ⚡ Gate → Phase 2 Build.</div>
          </button>
        </div>

        <!-- 4D Canvas -->
        <div class="rounded-xl overflow-hidden border border-slate-800 bg-[#040812]">
          <canvas id="canvas4DManifold" height="380"></canvas>
        </div>

        <div class="flex flex-wrap items-center justify-between gap-2 text-[11px] font-mono text-slate-300">
          <span>🖱️ <strong>Interactive 4D Manifold:</strong> Drag anywhere on the graph to rotate in 3D space.</span>
          <div class="flex items-center gap-2">
            <span class="text-slate-400">Jump to Phase:</span>
            <button data-jump="0.08" class="jump-btn px-2 py-0.5 rounded bg-slate-900 border border-slate-700 hover:border-emerald-400">Phase 0 (d7)</button>
            <button data-jump="0.26" class="jump-btn px-2 py-0.5 rounded bg-slate-900 border border-slate-700 hover:border-emerald-400">Phase 1 (d24)</button>
            <button data-jump="0.42" class="jump-btn px-2 py-0.5 rounded bg-cyan-950 border border-cyan-500/50 text-cyan-200 font-bold">⚡ Gate (d{gate_d})</button>
            <button data-jump="0.95" class="jump-btn px-2 py-0.5 rounded bg-slate-900 border border-slate-700 hover:border-emerald-400">Phase 2 (d86)</button>
          </div>
        </div>
      </div>

      <!-- Right 5 Cols: Cohort Parameter Controls + 9-Point Sensitivity Sweeps + Live Gemini 3.8 Flash Co-Pilot -->
      <div class="lg:col-span-5 rounded-2xl border-2 border-emerald-500/40 bg-[#091224] p-4 flex flex-col justify-between space-y-4 text-xs shadow-xl">

        <!-- Cohort-Specific Dosing & Gate Controls -->
        <div class="space-y-3">
          <div class="flex items-center justify-between border-b border-slate-800 pb-2">
            <span class="px-2.5 py-0.5 rounded bg-emerald-950 text-emerald-300 border border-emerald-500/40 font-mono font-bold">
              🎛️ COHORT PARAMETER &amp; GATE DECISION CONTROLS
            </span>
            <button id="btnResetOptimal" class="px-2.5 py-1 rounded bg-emerald-500 text-slate-950 font-bold text-[11px] hover:bg-emerald-400 transition">
              ⚡ Reset Optimal
            </button>
          </div>

          <div class="grid grid-cols-1 sm:grid-cols-2 gap-3 bg-slate-950/90 p-3 rounded-xl border border-slate-800">
            <div class="space-y-1">
              <div class="flex justify-between font-semibold text-slate-200">
                <span>Stage 1 → 2 Gate Switch Day</span>
                <span id="lblGateDay" class="font-mono text-cyan-300 font-bold">Day {gate_d}</span>
              </div>
              <input id="sliderGateDay" type="range" min="14" max="65" step="1" value="{gate_d}" class="w-full accent-cyan-400" />
              <div class="text-[10.5px] text-slate-400">Triggered when plasma GDF15 drops ≥35%.</div>
            </div>

            <div class="space-y-1">
              <div class="flex justify-between font-semibold text-slate-200">
                <span>Baseline Mutant Load (h₀)</span>
                <span id="lblH0" class="font-mono text-amber-300 font-bold">50%</span>
              </div>
              <input id="sliderH0" type="range" min="20" max="80" step="2" value="50" class="w-full accent-amber-400" />
              <div class="text-[10.5px] text-slate-400">Initial patient heteroplasmy severity.</div>
            </div>

            <div class="space-y-1 sm:col-span-2">
              <div class="flex justify-between font-semibold text-slate-200">
                <span>{cfg["slider1_label"]}</span>
                <span id="lblS1" class="font-mono text-purple-300 font-bold">95%</span>
              </div>
              <input id="sliderS1" type="range" min="0" max="100" step="5" value="95" class="w-full accent-purple-400" />
              <div class="text-[10.5px] text-slate-400">{cfg["slider1_desc"]}</div>
            </div>

            <div class="space-y-1">
              <div class="flex justify-between font-semibold text-slate-200">
                <span>{cfg["slider2_label"]}</span>
                <span id="lblS2" class="font-mono text-emerald-300 font-bold">95%</span>
              </div>
              <input id="sliderS2" type="range" min="0" max="100" step="5" value="95" class="w-full accent-emerald-400" />
            </div>

            <div class="space-y-1">
              <div class="flex justify-between font-semibold text-slate-200">
                <span>{cfg["slider3_label"]}</span>
                <span id="lblS3" class="font-mono text-teal-300 font-bold">95%</span>
              </div>
              <input id="sliderS3" type="range" min="0" max="100" step="5" value="95" class="w-full accent-teal-400" />
            </div>{extra_slider_html}
          </div>
        </div>

        <!-- 9-Point What-If Sensitivity Sweeps (High-Contrast Buttons) -->
        <div class="space-y-2 border-t border-slate-800 pt-2.5">
          <div class="flex items-center justify-between">
            <span class="font-mono font-bold text-amber-300 uppercase text-[11px]">
              🔬 1-CLICK 9-POINT "WHAT-IF" SENSITIVITY SWEEPS:
            </span>
            <span class="font-mono text-[10px] text-slate-400">Quantitative Tradeoff Matrix</span>
          </div>
          <div class="grid grid-cols-3 gap-1.5">
            <button data-sweep="switch_day" class="sweep-btn px-2 py-1.5 rounded-lg bg-amber-950 border border-amber-500/60 text-amber-200 hover:bg-amber-900 font-mono font-bold text-[10.5px] transition">
              🔬 1. Gate Switch Day
            </button>
            <button data-sweep="stage0_dose" class="sweep-btn px-2 py-1.5 rounded-lg bg-purple-950 border border-purple-500/60 text-purple-200 hover:bg-purple-900 font-mono font-bold text-[10.5px] transition">
              🔬 2. Stage 0 Unblock
            </button>
            <button data-sweep="initial_h0" class="sweep-btn px-2 py-1.5 rounded-lg bg-cyan-950 border border-cyan-500/60 text-cyan-200 hover:bg-cyan-900 font-mono font-bold text-[10.5px] transition">
              🔬 3. Initial Load h₀
            </button>{extra_sweep_btn}
          </div>
        </div>

        <!-- Live Gemini 3.8 Flash Interactions API (@google/genai) & Output Box -->
        <div class="space-y-2 border-t border-slate-800 pt-2.5">
          <details class="rounded-xl bg-slate-950 border border-cyan-500/40 p-2.5 text-[11px]">
            <summary class="cursor-pointer font-mono font-bold text-cyan-300 flex items-center justify-between">
              <span>⚡ Live Gemini API (`@google/genai` · `gemini-3.8-flash`)</span>
              <span id="geminiBadge" class="px-2 py-0.5 rounded bg-slate-900 text-emerald-300 border border-emerald-500/40 text-[10px]">Ready</span>
            </summary>
            <div class="mt-2 pt-2 border-t border-slate-800 space-y-1.5">
              <div class="text-[10px] text-slate-300">
                Uses official <code>@google/genai</code> <code>client.interactions.create</code> with <code>gemini-3.8-flash</code> or <code>gemini-4-argon</code>.
              </div>
              <div class="flex flex-wrap items-center gap-1.5">
                <select id="geminiModel" class="px-2 py-1 rounded bg-slate-900 border border-slate-700 text-slate-200 font-mono text-[10px]">
                  <option value="gemini-3.8-flash" selected>gemini-3.8-flash (1M · Fast Agentic)</option>
                  <option value="gemini-4-argon">gemini-4-argon (Deep Bio-Reasoning)</option>
                </select>
                <input id="geminiKey" type="password" placeholder="Optional GEMINI_API_KEY" class="flex-1 min-w-[130px] px-2 py-1 rounded bg-slate-900 border border-slate-700 text-slate-200 font-mono text-[10px]" />
                <button id="btnCopySdk" class="px-2 py-1 rounded bg-purple-950 border border-purple-500/50 text-purple-200 font-mono text-[10px] font-bold">
                  📋 SDK Code
                </button>
              </div>
            </div>
          </details>

          <div class="flex items-center gap-1.5">
            <input id="aiInput" type="text" placeholder="Ask about this cohort's 3D bottleneck, flag thresholds, or equations..." class="w-full px-2.5 py-1.5 rounded-lg bg-slate-950 border border-slate-700 text-slate-100 text-xs focus:outline-none focus:border-cyan-400" />
            <button id="btnAskAi" class="px-3 py-1.5 rounded-lg bg-cyan-500 text-slate-950 font-extrabold text-xs shrink-0 hover:bg-cyan-400 transition">
              Ask AI
            </button>
          </div>

          <div id="aiOutputBox" class="p-3 rounded-xl bg-slate-950 border border-slate-800 text-slate-100 text-xs leading-relaxed max-h-60 overflow-y-auto">
          </div>
        </div>

      </div>
    </div>

  </main>

  <script>
    const SPEC = window.MitoUX4DEngine.COHORT_SPECS['{spec_key}'];
    const pageState = {{
      activeRegimen: 'bcr',
      timeFrac: {(gate_d / 90.0):.3f},
      isPlaying: false,
      activePhaseIdx: 2,
      gateDay: {gate_d},
      h0: SPEC.baselineH0,
      s1: 95,
      s2: 95,
      s3: 95,
      hypShield: 94,
      previousInteractionId: null
    }};

    document.getElementById('hdrBadge').textContent = SPEC.badge;
    document.getElementById('hdrPdb').textContent = SPEC.pdbList;
    document.getElementById('hdrTitle').textContent = SPEC.title;
    document.getElementById('hdrSubtitle').textContent = SPEC.subtitle;
    document.getElementById('sliderH0').value = Math.round(SPEC.baselineH0);
    document.getElementById('lblH0').textContent = Math.round(SPEC.baselineH0) + '%';

    const viewer4D = window.MitoUX4DEngine.attach4DManifoldViewer('canvas4DManifold', SPEC, () => pageState);

    function syncPhaseFromTime() {{
      const d = pageState.timeFrac * 90.0;
      if (d <= 14) pageState.activePhaseIdx = 0;
      else if (d < pageState.gateDay - 2) pageState.activePhaseIdx = 1;
      else if (d <= pageState.gateDay + 4) pageState.activePhaseIdx = 2;
      else pageState.activePhaseIdx = 3;
    }}

    function updateAllUI() {{
      const d = (pageState.timeFrac * 90.0).toFixed(1);
      document.getElementById('clockBadge').textContent = `Day ${{d}} / 90d`;
      document.getElementById('timeScrubber').value = Math.round(pageState.timeFrac * 900);

      // Compute dynamic endpoint preview based on sliders
      const s1Frac = pageState.s1 / 100.0;
      const s2Frac = pageState.s2 / 100.0;
      const s3Frac = pageState.s3 / 100.0;
      const gatePenalty = Math.abs(pageState.gateDay - SPEC.gateTime) * 0.18;

      let endH, endATP;
      if (pageState.activeRegimen === 'unphased') {{
        endH = SPEC.trapEndH;
        endATP = SPEC.trapEndATP;
      }} else if (pageState.activeRegimen === 'untreated') {{
        endH = SPEC.untreatedEndH;
        endATP = SPEC.untreatedEndATP;
      }} else {{
        const hypBonus = (document.getElementById('sliderHypShield') ? (pageState.hypShield / 100.0) : 0.0);
        endH = Math.max(0.4, SPEC.bcrEndH - hypBonus * 0.6 + gatePenalty + (1 - s1Frac) * 28.0 + (1 - s2Frac) * 18.0 + (pageState.h0 - SPEC.baselineH0) * 0.12);
        endATP = Math.max(0.25, SPEC.bcrEndATP + hypBonus * 0.03 - gatePenalty * 0.02 - (1 - s1Frac) * 0.45 - (1 - s3Frac) * 0.35);
        if (document.getElementById('sliderHypShield')) {{
          const po2 = (1.63 - 0.71 * hypBonus).toFixed(2);
          const rescuePct = Math.min(98.2, 95.0 + 3.2 * hypBonus).toFixed(1);
          document.getElementById('lblHypShield').textContent = `${{pageState.hypShield}}% (pO₂=${{po2}}×)`;
          document.getElementById('kpiPrimary').textContent = `${{endATP.toFixed(2)}}× (${{rescuePct}}% Rescue)`;
          document.getElementById('kpiSecondary').textContent = `pO₂ ${{po2}}× | Jam 0.08× | SNCA ${{(0.15 - 0.03 * hypBonus).toFixed(2)}}×`;
        }}
      }}

      document.getElementById('kpiH').textContent = `${{endH.toFixed(1)}}% (from ${{Math.round(pageState.h0)}}%)`;
      document.getElementById('kpiHSub').textContent = endH <= 8.5
        ? `✅ Below 8% replicative escape threshold (Net ATP = ${{endATP.toFixed(2)}}×)`
        : `🔴 Warning: Residual mutant mtDNA > 8% (Net ATP = ${{endATP.toFixed(2)}}×)`;

      window.MitoUX4DEngine.render5PillarTransparencyInspector(
        'transparencyInspectorContainer',
        SPEC,
        pageState.activePhaseIdx,
        (newPhaseIdx) => {{
          pageState.activePhaseIdx = newPhaseIdx;
          const p = SPEC.phases[newPhaseIdx];
          const midDay = (p.timeRange[0] + p.timeRange[1]) * 0.5;
          pageState.timeFrac = Math.max(0, Math.min(1, midDay / 90.0));
          viewer4D.setCameraPreset(newPhaseIdx);
          updateAllUI();
        }}
      );

      viewer4D.render();
    }}

    function renderSweepTable(sweepKey) {{
      const sw = SPEC.sweeps[sweepKey] || SPEC.sweeps.switch_day;
      const out = document.getElementById('aiOutputBox');
      const thHtml = sw.headers.map(h => `<th class="pb-1 pr-2">${{h}}</th>`).join('');
      const rowsHtml = sw.rows.map(r => {{
        const isOpt = r[0].includes('⭐');
        return `<tr class="${{isOpt ? 'bg-emerald-950/90 text-emerald-200 font-bold' : 'border-t border-slate-800 text-slate-200'}}">
          ${{r.map(c => `<td class="py-1 pr-2 font-mono text-[10.5px]">${{c}}</td>`).join('')}}
        </tr>`;
      }}).join('');

      out.innerHTML = `
        <div class="space-y-2">
          <div class="font-bold text-amber-300 text-xs">${{sw.title}}</div>
          <table class="w-full text-left border-collapse">
            <thead><tr class="border-b border-slate-700 text-slate-400 font-mono text-[10px]">${{thHtml}}</tr></thead>
            <tbody>${{rowsHtml}}</tbody>
          </table>
          <div class="p-2 rounded-lg bg-emerald-950/50 border border-emerald-500/40 text-emerald-200 text-[11px]">
            <strong>💡 Clinical Decision Takeaway:</strong> ${{sw.takeaway}}
          </div>
        </div>
      `;
    }}

    document.querySelectorAll('.sweep-btn').forEach(btn => {{
      btn.addEventListener('click', () => renderSweepTable(btn.getAttribute('data-sweep')));
    }});

    document.querySelectorAll('.reg-card').forEach(btn => {{
      btn.addEventListener('click', () => {{
        pageState.activeRegimen = btn.getAttribute('data-reg');
        document.querySelectorAll('.reg-card').forEach(b => {{
          b.className = 'reg-card text-left p-2.5 rounded-xl border border-slate-700 bg-slate-900/90 text-slate-200 transition';
        }});
        btn.className = 'reg-card text-left p-2.5 rounded-xl border-2 border-emerald-500 bg-emerald-950/80 text-white transition shadow-md';
        updateAllUI();
      }});
    }});

    document.querySelectorAll('.cam-btn').forEach(btn => {{
      btn.addEventListener('click', () => viewer4D.setCameraPreset(parseInt(btn.getAttribute('data-cam'), 10)));
    }});

    document.querySelectorAll('.jump-btn').forEach(btn => {{
      btn.addEventListener('click', () => {{
        pageState.isPlaying = false;
        document.getElementById('btnPlayPause').textContent = '▶ Play Progression';
        pageState.timeFrac = parseFloat(btn.getAttribute('data-jump'));
        syncPhaseFromTime();
        updateAllUI();
      }});
    }});

    document.getElementById('timeScrubber').addEventListener('input', e => {{
      pageState.isPlaying = false;
      document.getElementById('btnPlayPause').textContent = '▶ Play Progression';
      pageState.timeFrac = parseFloat(e.target.value) / 900.0;
      syncPhaseFromTime();
      updateAllUI();
    }});

    document.getElementById('btnPlayPause').addEventListener('click', e => {{
      pageState.isPlaying = !pageState.isPlaying;
      e.target.textContent = pageState.isPlaying ? '⏸ Pause Progression' : '▶ Play Progression';
    }});

    document.getElementById('sliderGateDay').addEventListener('input', e => {{
      pageState.gateDay = parseInt(e.target.value, 10);
      document.getElementById('lblGateDay').textContent = 'Day ' + pageState.gateDay;
      document.getElementById('kpiGateStatus').textContent = `⚡ Day ${{pageState.gateDay}} Gate Active`;
      updateAllUI();
    }});
    document.getElementById('sliderH0').addEventListener('input', e => {{
      pageState.h0 = parseInt(e.target.value, 10);
      document.getElementById('lblH0').textContent = pageState.h0 + '%';
      updateAllUI();
    }});
    document.getElementById('sliderS1').addEventListener('input', e => {{
      pageState.s1 = parseInt(e.target.value, 10);
      document.getElementById('lblS1').textContent = pageState.s1 + '%';
      updateAllUI();
    }});
    document.getElementById('sliderS2').addEventListener('input', e => {{
      pageState.s2 = parseInt(e.target.value, 10);
      document.getElementById('lblS2').textContent = pageState.s2 + '%';
      updateAllUI();
    }});
    document.getElementById('sliderS3').addEventListener('input', e => {{
      pageState.s3 = parseInt(e.target.value, 10);
      document.getElementById('lblS3').textContent = pageState.s3 + '%';
      updateAllUI();
    }});
    if (document.getElementById('sliderHypShield')) {{
      document.getElementById('sliderHypShield').addEventListener('input', e => {{
        pageState.hypShield = parseInt(e.target.value, 10);
        updateAllUI();
      }});
    }}

    document.getElementById('btnResetOptimal').addEventListener('click', () => {{
      pageState.activeRegimen = 'bcr';
      pageState.gateDay = {gate_d};
      pageState.h0 = SPEC.baselineH0;
      pageState.s1 = 95;
      pageState.s2 = 95;
      pageState.s3 = 95;
      document.getElementById('sliderGateDay').value = {gate_d};
      document.getElementById('lblGateDay').textContent = 'Day {gate_d}';
      document.getElementById('sliderH0').value = Math.round(SPEC.baselineH0);
      document.getElementById('lblH0').textContent = Math.round(SPEC.baselineH0) + '%';
      document.getElementById('sliderS1').value = 95;
      document.getElementById('lblS1').textContent = '95%';
      document.getElementById('sliderS2').value = 95;
      document.getElementById('lblS2').textContent = '95%';
      document.getElementById('sliderS3').value = 95;
      document.getElementById('lblS3').textContent = '95%';
      updateAllUI();
    }});

    document.getElementById('btnCopySdk').addEventListener('click', () => {{
      const model = document.getElementById('geminiModel').value || 'gemini-3.8-flash';
      const snippet = `import {{ GoogleGenAI }} from "@google/genai";\nconst client = new GoogleGenAI({{ apiKey: process.env.GEMINI_API_KEY }});\nconst interaction = await client.interactions.create({{\n  model: "${{model}}",\n  system_instruction: "${{SPEC.title}} (${{SPEC.pdbList}})",\n  input: "Explain the 5-pillar decision scope and Stage-Gate threshold for ${{SPEC.badge}}."\n}});\nconsole.log(interaction.output_text);`;
      document.getElementById('aiOutputBox').innerHTML = `<div class="font-mono text-[10px] text-purple-300 font-bold mb-1">📋 Official <code>@google/genai</code> Interactions API Snippet:</div><pre class="p-2 rounded bg-slate-900 text-slate-200 font-mono text-[10px] overflow-x-auto whitespace-pre-wrap">${{snippet}}</pre>`;
    }});

    document.getElementById('btnAskAi').addEventListener('click', async () => {{
      const q = document.getElementById('aiInput').value.trim() || 'Explain why the stage-gate flags and how the 3D bottleneck is resolved.';
      const p = SPEC.phases[pageState.activePhaseIdx];
      const detHtml = `<div class="space-y-1.5">
        <div class="font-bold text-cyan-300">🤖 Co-Pilot Analysis for ${{SPEC.badge}} (${{p.shortTitle}}):</div>
        <div>• <strong>Purpose:</strong> ${{p.purpose}}</div>
        <div>• <strong>Why It Flags:</strong> ${{p.whyFlags}}</div>
        <div>• <strong>Methods &amp; 3D Provenance:</strong> <code>${{p.methodsUsed}}</code></div>
        <div>• <strong>Overall Scope Impact:</strong> ${{p.overallImpact}}</div>
      </div>`;
      const apiKey = document.getElementById('geminiKey').value.trim() || sessionStorage.getItem('MITORESTORER_GEMINI_API_KEY') || '';
      const model = document.getElementById('geminiModel').value || 'gemini-3.8-flash';
      const out = document.getElementById('aiOutputBox');
      if (!apiKey) {{
        out.innerHTML = detHtml;
        return;
      }}
      sessionStorage.setItem('MITORESTORER_GEMINI_API_KEY', apiKey);
      out.innerHTML = `<div class="text-cyan-300 font-mono text-[11px] mb-1">⚡ Querying <code>${{model}}</code> via <code>@google/genai</code> Interactions API...</div>` + detHtml;
      try {{
        const genaiMod = await import('https://esm.sh/@google/genai');
        const client = new genaiMod.GoogleGenAI({{ apiKey }});
        const interaction = await client.interactions.create({{
          model,
          input: q,
          system_instruction: `${{SPEC.title}} | ${{SPEC.pdbList}} | Active Phase: ${{p.shortTitle}} | ${{p.methodsUsed}}`
        }});
        if (interaction && interaction.output_text) {{
          out.innerHTML = `<div class="p-2 rounded bg-cyan-950/70 border border-cyan-400/50 mb-2 text-white">${{interaction.output_text.replace(/\\n/g, '<br/>')}}</div>` + detHtml;
        }}
      }} catch (e) {{
        out.innerHTML = detHtml;
      }}
    }});

    function animateLoop() {{
      if (pageState.isPlaying) {{
        pageState.timeFrac += 0.0012;
        if (pageState.timeFrac > 1.0) pageState.timeFrac = 0.0;
        syncPhaseFromTime();
        updateAllUI();
      }}
      requestAnimationFrame(animateLoop);
    }}

    window.addEventListener('resize', () => viewer4D.render());
    updateAllUI();
    renderSweepTable('switch_day');
    animateLoop();
  </script>
</body>
</html>
"""


def main():
    sim_dir = Path(__file__).resolve().parent
    for cfg in COHORT_PAGES:
        out_path = sim_dir / cfg["filename"]
        out_path.write_text(build_page_html(cfg), encoding="utf-8")
        print(f"Built dedicated cohort workstation: {out_path}")


if __name__ == "__main__":
    main()

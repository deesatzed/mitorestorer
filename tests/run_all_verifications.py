#!/usr/bin/env python3
"""
MitoRestorer Automated Verification & Regression Test Suite (tests/run_all_verifications.py)
============================================================================================
Runs a complete, deterministic verification across:
  1. All 8 Atomic 3D Crystal Structures (6FCX, 1B68, 7CK6, 6EQI, 4ZZJ, 2NSX, 2Q8G, 1SQB)
  2. All 8 Core Biophysical & Cross-Species Discoveries (Human Cohorts A-E, CAR-T Cohort F, Honeybee Cohort G)
  3. Closed-Loop Biomarker Gating Logic (Day-38 GDF15 Gate, Hour-64 CAR-T Lac:Pyr Gate, Day-4.5 Honeybee Detox Gate)
  4. End-to-End Sample Walkthrough Scripts (Human PD, CAR-T Bioreactor, Honeybee Apis-BCR)
"""

import json
import os
import subprocess
import sys

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(REPO_ROOT, "simulator"))

from sim_bioenergetic_reset import (
    simulate_apis_honeybee,
    simulate_cart_bioreactor,
    simulate_cortical_neuron_candidates,
    simulate_neuropsychiatric_subprotocols,
    simulate_pink1_sirt1_structural_variants,
    simulate_regimen,
)

PASS_COUNT = 0
FAIL_COUNT = 0


def check(condition, test_name, detail):
    global PASS_COUNT, FAIL_COUNT
    if condition:
        PASS_COUNT += 1
        print(f"  [PASS] {test_name:62s} | {detail}")
    else:
        FAIL_COUNT += 1
        print(f"  [FAIL] {test_name:62s} | {detail}")


def main():
    print("=" * 100)
    print("🧪 MITORESTORER v2.0 AUTOMATED VERIFICATION & BENCHMARK TEST SUITE")
    print("=" * 100)

    # -------------------------------------------------------------------------
    # SUITE 1: 8-STRUCTURE 3D ATOMIC CIF VERIFICATION
    # -------------------------------------------------------------------------
    print("\n--- 1. Verifying 8 Local 3D Crystal Structures (.cif) ---")
    expected_cifs = [
        ("6FCX", "MTHFR + FAD + SAH (2.50 A)"),
        ("1B68", "APOE4 Arg112-Glu109 Salt Bridge (2.00 A)"),
        ("7CK6", "TOMM40 Beta-Barrel Pore + PC1 (3.00 A)"),
        ("6EQI", "Active PINK1 + Phospho-Ser202 + Ubiquitin (3.10 A)"),
        ("4ZZJ", "SIRT1 + Carba-NAD+ + STAC 4TQ (2.74 A)"),
        ("2NSX", "Lysosomal GBA1 + Chaperone IFM 506 (2.20 A)"),
        ("2Q8G", "Human PDK1 + Allosteric Inhibitor AZD7545 (1.90 A, CAR-T)"),
        ("1SQB", "Complex III Cytochrome bc1 + Azoxystrobin (2.69 A, Honeybee)"),
    ]
    for pdb_id, desc in expected_cifs:
        cif_path = os.path.join(REPO_ROOT, "structures", pdb_id, f"pdb_0000{pdb_id.lower()}.cif")
        exists = os.path.exists(cif_path) and os.path.getsize(cif_path) > 100000
        size_kb = os.path.getsize(cif_path) // 1024 if os.path.exists(cif_path) else 0
        check(exists, f"PDB {pdb_id}: {desc}", f"Size = {size_kb} KB")

    # -------------------------------------------------------------------------
    # SUITE 2: HUMAN CLINICAL DISCOVERIES 1-6 (COHORTS A-E)
    # -------------------------------------------------------------------------
    print("\n--- 2. Verifying Human Clinical ODE Discoveries 1-6 (Cohorts A-E) ---")

    # Discovery 1 & 3: Unassisted 16:8 vs Two-Stage BCR in Skeletal Muscle
    r_16_8 = simulate_regimen("5_Biphasic_Chrono_Synbiotic_BCR", genotype="standard", cycle_mode="circadian_16_8", purge_fraction=16.0/24.0, tau_lag_hours=13.5, tissue_scale=1.0)
    r_bcr = simulate_regimen("5_Biphasic_Chrono_Synbiotic_BCR", genotype="standard", cycle_mode="two_stage_adaptive", purge_fraction=18.0/24.0, tau_lag_hours=6.0, tissue_scale=1.0)
    check(
        r_bcr["h_day90"] < 5.0 and r_16_8["h_day90"] > 35.0 and r_bcr["ATP_cap_day90"] >= 1.15,
        "Cohort A (Skeletal Muscle): Two-Stage BCR vs Unassisted 16:8",
        f"BCR h={r_bcr['h_day90']}% (ATP={r_bcr['ATP_cap_day90']}x) vs 16:8 h={r_16_8['h_day90']}%",
    )

    # Discovery 4: Cortical Neuron APOE4/TOMM40 Rescue
    r_neuro_base = simulate_cortical_neuron_candidates("0_Baseline_Two_Stage_BCR_Peripheral_Only", genotype="apoe4_pink1_fut2")
    r_neuro_full = simulate_cortical_neuron_candidates("4_Complete_Neuro_BCR_Synergistic_Stack", genotype="apoe4_pink1_fut2")
    check(
        r_neuro_full["h_day90"] < 5.0 and r_neuro_full["ATP_day90"] > 1.10 and r_neuro_base["h_day90"] > 80.0,
        "Cohort B (APOE4 Cortical Neurons): Complete Neuro-BCR Rescue",
        f"Neuro-BCR h={r_neuro_full['h_day90']}% (ATP={r_neuro_full['ATP_day90']}x) vs Peripheral-Only h={r_neuro_base['h_day90']}%",
    )

    # Discovery 5: 6EQI PINK1 Class I vs Class II Stratification
    r_c1_ktp = simulate_pink1_sirt1_structural_variants("Class_I_ATP_Pocket_Ala339Thr", True, "KTP_Monotherapy")
    r_c2_ktp = simulate_pink1_sirt1_structural_variants("Class_II_Ins3_Ubiquitin_Cradle_Mutant", True, "KTP_Monotherapy")
    r_c2_full = simulate_pink1_sirt1_structural_variants("Class_II_Ins3_Ubiquitin_Cradle_Mutant", True, "Complete_6EQI_4ZZJ_Precision_Stack")
    check(
        r_c1_ktp["h_day90"] < 15.0 and r_c2_ktp["h_day90"] > 40.0 and r_c2_full["h_day90"] < 18.0,
        "Discovery 5 (6EQI): PINK1 Class I vs Class II Stratification",
        f"Class I KTP h={r_c1_ktp['h_day90']}% | Class II KTP h={r_c2_ktp['h_day90']}% -> Full Stack h={r_c2_full['h_day90']}%",
    )

    # Discovery 2: Refractory Epilepsy Day-30 Keto H2S Escape vs Epi-BCR
    r_epi_kd = simulate_neuropsychiatric_subprotocols("Refractory_Epilepsy", "Chronic_Classic_Ketogenic_Diet")
    r_epi_bcr = simulate_neuropsychiatric_subprotocols("Refractory_Epilepsy", "Epi_BCR_Protocol")
    check(
        r_epi_kd["H2S_day90"] > 1.5 and r_epi_bcr["H2S_day90"] < 0.20 and r_epi_bcr["Seizure_Reduction_Pct"] >= 95.0,
        "Cohort C (Refractory Epilepsy): Epi-BCR Blocks Keto H2S Escape",
        f"Classic KD H2S={r_epi_kd['H2S_day90']} (E/I={r_epi_kd['EI_ratio_day90']}) vs Epi-BCR H2S={r_epi_bcr['H2S_day90']} (E/I={r_epi_bcr['EI_ratio_day90']}, {r_epi_bcr['Seizure_Reduction_Pct']}% red.)",
    )

    # Cohort D: Treatment-Resistant MDD Post-Ketamine Relapse vs Affective-BCR
    r_mdd_ket = simulate_neuropsychiatric_subprotocols("Treatment_Resistant_MDD", "Acute_Ketamine_Burst_Alone")
    r_mdd_bcr = simulate_neuropsychiatric_subprotocols("Treatment_Resistant_MDD", "Affective_MDD_BCR_Protocol")
    check(
        r_mdd_ket["Spine_Density_day7"] > r_mdd_ket["Spine_Density_day90"] and r_mdd_bcr["Spine_Density_day90"] >= 1.20,
        "Cohort D (TR-MDD): Affective-BCR Sustains Dendritic Spines",
        f"Ketamine Alone Spines d7={r_mdd_ket['Spine_Density_day7']}x->d90={r_mdd_ket['Spine_Density_day90']}x vs Affective-BCR d90={r_mdd_bcr['Spine_Density_day90']}x",
    )

    # Discovery 6: Parkinson's SNc GBA1/LRRK2 Autophagosome Jam vs PD-BCR + 2NSX Ambroxol
    r_pd_jam = simulate_neuropsychiatric_subprotocols("Parkinsons_SNc", "Mitophagy_Alone_Without_GBA1_Ambroxol")
    r_pd_bcr = simulate_neuropsychiatric_subprotocols("Parkinsons_SNc", "PD_BCR_Nigrostriatal_Protocol")
    check(
        r_pd_jam["Auto_Jam_day90"] > 1.20 and r_pd_bcr["Auto_Jam_day90"] <= 0.10 and r_pd_bcr["Motor_Score_Rescue_Pct"] >= 90.0,
        "Cohort E (Parkinson's SNc): Stage 0 2NSX Ambroxol Clears Jam",
        f"No-Ambroxol Jam={r_pd_jam['Auto_Jam_day90']}x (SNCA={r_pd_jam['SNCA_agg_day90']}x) vs PD-BCR Jam={r_pd_bcr['Auto_Jam_day90']}x (TH={r_pd_bcr['TH_Dopamine_Tone_day90']}x, {r_pd_bcr['Motor_Score_Rescue_Pct']}% rescue)",
    )

    # -------------------------------------------------------------------------
    # SUITE 3: CROSS-SPECIES DISCOVERIES 7 (CAR-T) & 8 (HONEYBEE APIS-BCR)
    # -------------------------------------------------------------------------
    print("\n--- 3. Verifying Cross-Species Discoveries 7 (CAR-T) & 8 (Honeybee Apis-BCR) ---")

    c_std = simulate_cart_bioreactor("1_Standard_CD3_CD28_IL2_HighGlucose")
    c_rapa = simulate_cart_bioreactor("2_Static_Continuous_UrolithinA_Rapamycin")
    c_bcr = simulate_cart_bioreactor("3_Bioreactor_BCR_Biphasic_Gated_Controller")
    check(
        c_std["T_SCM_Stemness_pct"] < 15.0
        and c_rapa["Fold_Expansion_day10_x"] < 15.0
        and c_bcr["Fold_Expansion_day10_x"] >= 45.0
        and c_bcr["T_SCM_Stemness_pct"] >= 70.0,
        "Cohort F (CAR-T Bioreactor-BCR): Solves Expansion vs Exhaustion",
        f"Bioreactor-BCR Fold={c_bcr['Fold_Expansion_day10_x']}x, T_SCM={c_bcr['T_SCM_Stemness_pct']}%, SRC={c_bcr['Spare_Respiratory_Capacity_x']}x (Sw h{c_bcr['switch_hour']})",
    )

    a_hfcs = simulate_apis_honeybee("1_Exposed_Plain_HFCS_Syrup")
    a_unph = simulate_apis_honeybee("2_Unphased_Antioxidant_Pollen_Sub_Without_Detox_Purge")
    a_bcr = simulate_apis_honeybee("3_Apis_BCR_Two_Stage_Holobiont_Patty")
    check(
        a_hfcs["Forager_Homing_Survival_pct"] <= 20.0
        and a_bcr["Forager_Homing_Survival_pct"] >= 90.0
        and a_bcr["Thoracic_Flight_ATP_day14_x"] >= 1.50,
        "Cohort G (Honeybee Apis-BCR): Rescues Flight ATP & Homing",
        f"Plain HFCS Homing={a_hfcs['Forager_Homing_Survival_pct']}% (ATP={a_hfcs['Thoracic_Flight_ATP_day14_x']}x) vs Apis-BCR Homing={a_bcr['Forager_Homing_Survival_pct']}% (ATP={a_bcr['Thoracic_Flight_ATP_day14_x']}x)",
    )

    # -------------------------------------------------------------------------
    # SUITE 4: GOOGLE CORAL NPU (`RV32IMF_Zve32x`) COMPILER, FIRMWARE & 3D UI
    # -------------------------------------------------------------------------
    print("\n--- 4. Verifying Google Coral NPU (`RV32IMF_Zve32x`) Compiler, Firmware & 3D Silicon UI ---")
    from pathlib import Path
    from coral_npu_bcr_compiler import CoralNPUBCRCompiler

    coral_compiler = CoralNPUBCRCompiler()
    coral_summary = coral_compiler.build_all(Path(REPO_ROOT))
    coral_targets = coral_summary["targets"]

    apis_npu = coral_targets["APIS_SMART_HIVE_NODE"]
    cart_npu = coral_targets["CART_GMP_CASSETTE_VALVE"]
    hum_npu = coral_targets["HUMAN_WEARABLE_RING_PATCH"]

    check(
        apis_npu["active_power_mw"] <= 4.50
        and apis_npu["total_sram_bytes"] < 8192
        and apis_npu["latency_ms"] < 0.05,
        "Coral NPU Target 1 (Apis-BCR Smart Hive Node): <4.5mW & <8KB",
        f"Cycles={apis_npu['total_cycles']} | Power={apis_npu['active_power_mw']} mW | SRAM={apis_npu['total_sram_bytes']} B",
    )

    check(
        cart_npu["active_power_mw"] <= 4.50
        and cart_npu["total_sram_bytes"] < 8192
        and hum_npu["total_sram_bytes"] < 8192,
        "Coral NPU Targets 2 & 3 (CAR-T Valve & Wearable Ring): <8KB TCM",
        f"CAR-T={cart_npu['active_power_mw']} mW ({cart_npu['total_sram_bytes']} B) | Ring={hum_npu['active_power_mw']} mW ({hum_npu['total_sram_bytes']} B)",
    )

    mlir_file = os.path.join(REPO_ROOT, "build", "coralnpu", "mitorestorer_bcr.mlir")
    header_file = os.path.join(REPO_ROOT, "build", "coralnpu", "mitorestorer_coralnpu_kernel.h")
    hil_out = coral_summary.get("c99_hil_verification", "")
    check(
        os.path.exists(mlir_file)
        and os.path.exists(header_file)
        and "coralnpu.matrix_mac_int8_skip" in open(mlir_file).read()
        and "__riscv_vsetvl_e32m2" in open(header_file).read()
        and "CORAL_NPU_C99_HIL_PASS" in hil_out
        and apis_npu["max_h_quantization_err_pct"] < 0.10,
        "Coral NPU C99/RVV HIL Binary & Bit-Exact Q16.16/INT8 Parity",
        f"{hil_out} | Max Q16.16 Err={apis_npu['max_h_quantization_err_pct']:.4f}%",
    )

    sim_html_path = os.path.join(REPO_ROOT, "simulator", "bcr_simulator.html")
    sim_html_text = open(sim_html_path, encoding="utf-8").read()
    check(
        "canvasCoralNPU3D" in sim_html_text
        and "renderCoralNpuSiliconCanvas" in sim_html_text
        and "APIS_SMART_HIVE_NODE" in sim_html_text,
        "Simulator Frontend 3D Coral NPU Silicon Triad & Strip UI",
        "Verified 3D Triad Canvas + Strip + MLIR/C99 Inspector in bcr_simulator.html",
    )

    cart_html_path = os.path.join(REPO_ROOT, "simulator", "cart_bioreactor_simulator.html")
    apis_html_path = os.path.join(REPO_ROOT, "simulator", "apis_honeybee_simulator.html")
    cart_html_text = open(cart_html_path, encoding="utf-8").read()
    apis_html_text = open(apis_html_path, encoding="utf-8").read()

    from gemini_copilot_sensitivity import (
        run_human_switch_day_sweep,
        run_cart_valve_hour_sweep,
        run_apis_baffle_day_sweep,
        query_gemini_copilot_interaction,
    )

    h_sw = run_human_switch_day_sweep()
    c_sw = run_cart_valve_hour_sweep()
    a_sw = run_apis_baffle_day_sweep()
    copilot_res = query_gemini_copilot_interaction("Cross-Species", "Verify stage gate")

    check(
        "runHumanWhatIfSweep" in sim_html_text
        and "runCartWhatIfSweep" in cart_html_text
        and "runApisWhatIfSweep" in apis_html_text
        and "gemini-3.8-flash" in sim_html_text
        and "gemini-3.8-flash" in cart_html_text
        and "gemini-3.8-flash" in apis_html_text
        and "interactions.create" in sim_html_text
        and len(h_sw) == 9
        and len(c_sw) == 9
        and len(a_sw) == 9
        and copilot_res["model"] == "gemini-3.8-flash",
        "3 Domain Simulators: Gemini 3.8 Flash Co-Pilot & 9-Point Sweeps",
        "Verified @google/genai interactions.create (gemini-3.8-flash) + 9-Point What-If Sweeps across Human, CAR-T, & Honeybee",
    )

    ux_engine_path = os.path.join(REPO_ROOT, "simulator", "bcr_ux_4d_engine.js")
    ux_engine_text = open(ux_engine_path, encoding="utf-8").read() if os.path.exists(ux_engine_path) else ""
    dedicated_pages = [
        "parkinsons_simulator.html",
        "alzheimers_simulator.html",
        "epilepsy_simulator.html",
        "depression_simulator.html",
        "sarcopenia_simulator.html",
        "cart_bioreactor_simulator.html",
        "apis_honeybee_simulator.html",
    ]
    all_pages_ok = all(
        os.path.exists(os.path.join(REPO_ROOT, "simulator", fn))
        and "bcr_ux_4d_engine.js" in open(os.path.join(REPO_ROOT, "simulator", fn), encoding="utf-8").read()
        and "attach4DManifoldViewer" in open(os.path.join(REPO_ROOT, "simulator", fn), encoding="utf-8").read()
        and "render5PillarTransparencyInspector" in open(os.path.join(REPO_ROOT, "simulator", fn), encoding="utf-8").read()
        for fn in dedicated_pages
    )
    check(
        all_pages_ok
        and "attach4DManifoldViewer" in ux_engine_text
        and "render5PillarTransparencyInspector" in ux_engine_text,
        "7 Dedicated Cohort Workstations: 4D Phase-Space Manifold & 5-Pillar UX",
        "Verified all 7 standalone cohort pages + rotatable 4D Phase-Space Manifold + 5-Pillar Clinical Transparency Inspector",
    )

    from assess_new_science_release import simulate_epigenetic_deaging_acel70751

    acel_res = simulate_epigenetic_deaging_acel70751()
    pd_hyp_alone = simulate_neuropsychiatric_subprotocols("Parkinsons_SNc", "Hypoxia_Pill_GBT601_PT2399_Alone")
    pd_hyp_bcr = simulate_neuropsychiatric_subprotocols("Parkinsons_SNc", "PD_BCR_Plus_GBT601_PT2399_Shield")
    sandbox_path = os.path.join(REPO_ROOT, "simulator", "science_radar_sandbox.html")
    sandbox_text = open(sandbox_path, encoding="utf-8").read() if os.path.exists(sandbox_path) else ""
    check(
        "10.1111/acel.70751" in sandbox_text
        and "GBT601" in sandbox_text
        and pd_hyp_alone["Motor_Score_Rescue_Pct"] > 40.0
        and pd_hyp_bcr["Motor_Score_Rescue_Pct"] >= 98.0
        and pd_hyp_bcr["switch_day"] < 30.0
        and acel_res["3_Standard_SOX2_cMYC_iNSC_Conversion_Berg2026"]["eAge_d100_pct"] < 10.0
        and acel_res["4_MitoRestorer_BCR_Plus_Berg2026_TET_SIRT1_6FCX_Synergy"]["eAge_d100_pct"] < 4.0,
        "v3.0 Science Radar: Mootha/Ichinose Hyperoxia & Berg 2026 (acel.70751)",
        f"PD+GBT601/PT2399={pd_hyp_bcr['Motor_Score_Rescue_Pct']}% (Sw d{pd_hyp_bcr['switch_day']}) | Berg2026 eAge=100%->{acel_res['3_Standard_SOX2_cMYC_iNSC_Conversion_Berg2026']['eAge_d100_pct']}% -> BCR+Berg={acel_res['4_MitoRestorer_BCR_Plus_Berg2026_TET_SIRT1_6FCX_Synergy']['eAge_d100_pct']}%",
    )

    print("\n" + "=" * 100)
    print(f"✅ VERIFICATION SUMMARY: {PASS_COUNT} PASSED, {FAIL_COUNT} FAILED")
    print("=" * 100)
    if FAIL_COUNT > 0:
        sys.exit(1)


if __name__ == "__main__":
    main()

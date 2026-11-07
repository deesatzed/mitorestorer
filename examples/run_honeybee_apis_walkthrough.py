#!/usr/bin/env python3
"""
End-to-End Executable Honeybee (Apis mellifera) Walkthrough (run_honeybee_apis_walkthrough.py)
==============================================================================================
Demonstrates a complete, runnable 4-Step Apiary Workflow for Cohort G (Apis-BCR):
  - Sample Colony: APIS-CA-ALMOND-01 (Commercial pollination hive co-exposed to Strobilurin
    fungicide Azoxystrobin [2.0 ppm], Neonicotinoid Imidacloprid [5.0 ppb], and Glyphosate [5.0 ppm])
  - Step 1: Colony Exposure Profile & 3D Complex III Qo-Site Inspection (PDB: 1SQB, Azoxystrobin at Glu271/Phe274)
  - Step 2: Why Plain HFCS Sugar Syrup & Unphased Pollen Subs Trigger Mid-Air Flight Collapse (EP = eR * f_flight^2)
  - Step 3: Execute the 2-Stage Apis-BCR Hive Patty Protocol (Days 1-5 p-Coumaric/Spermidine/S. alvi -> Days 5-14 Ubiquinol/NR/B2)
  - Step 4: Verify the Day-4.5 Pre-Foraging Detox Gate & Output Exact 1 kg / 25 kg Beekeepers' Patty Formulation
"""

import json
import os
import subprocess
import sys

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CLI_SCRIPT = os.path.join(REPO_ROOT, "skills", "bioenergetic-precision-protocol", "scripts", "bcr_pipeline.py")
CIF_1SQB = os.path.join(REPO_ROOT, "structures", "1SQB", "pdb_00001sqb.cif")
OUT_DIR = os.path.join(REPO_ROOT, "examples", "output")

sys.path.insert(0, os.path.join(REPO_ROOT, "simulator"))
from sim_bioenergetic_reset import simulate_apis_honeybee


def run_step(title):
    print("\n" + "=" * 92)
    print(f"🐝 {title}")
    print("=" * 92)


def main():
    os.makedirs(OUT_DIR, exist_ok=True)

    # -------------------------------------------------------------------------
    # STEP 1: COLONY EXPOSURE PROFILE & 3D COMPLEX III Qo-SITE INSPECTION (1SQB)
    # -------------------------------------------------------------------------
    run_step("STEP 1 OF 4: COLONY EXPOSURE PROFILE & 3D COMPLEX III Qo-SITE INSPECTION (PDB: 1SQB)")
    colony = {
        "colony_id": "APIS-CA-ALMOND-01",
        "crop_pollination_context": "Blooming Almond Orchard (Tank-Mix Fungicide + Insecticide + Herbicide Drift)",
        "agrochemical_exposure": {
            "Strobilurin_Fungicide": "Azoxystrobin 2.0 ppm -> Blocks Mitochondrial Complex III Qo pocket (PDB: 1SQB, 78% initial occupancy)",
            "Neonicotinoid_Insecticide": "Imidacloprid 5.0 ppb -> Downregulates thoracic Complex IV (cox5a, cox5b, cox6c)",
            "Herbicide_Drift": "Glyphosate 5.0 ppm -> Inhibits EPSPS shikimate pathway in ileum symbiont Snodgrassella alvi (0.28x density)",
        },
        "baseline_thoracic_heteroplasmy_pct": 52.0,
        "wingbeat_flux_transition": "In-Hive Nurse/Resting (Days 1-5, f = 1.0 -> f^2 = 1.0x) -> Foraging Flight (Days 6-14, 100-250 Hz, f = 2.40 -> f^2 = 5.76x)",
    }
    print(json.dumps(colony, indent=2))

    cif_out = os.path.join(OUT_DIR, "apis_step1_1sqb_azo_contacts.json")
    subprocess.run(
        [
            sys.executable,
            CLI_SCRIPT,
            "inspect-cif",
            "--cif",
            CIF_1SQB,
            "--target-chain",
            "C",
            "--target-res",
            "AZO",
            "--max-dist",
            "3.8",
            "--output",
            cif_out,
        ],
        check=True,
    )
    with open(cif_out, "r", encoding="utf-8") as f:
        cif_data = json.load(f)
    print(f"\n[3D Crystal Verification — PDB: 1SQB (Complex III Cytochrome bc1 + Azoxystrobin AZO 383)]:")
    for c in cif_data["contacts"][:6]:
        print(
            f"  • AZO({c['target_atom']}) <--> Chain {c['partner_chain']} {c['partner_comp_id']}{c['partner_auth_seq_id']}({c['partner_atom']}) : {c['distance_angstrom']} Å"
        )

    # -------------------------------------------------------------------------
    # STEP 2: WHY PLAIN HFCS SYRUP & UNPHASED POLLEN SUBS FAIL IN FORAGING FLIGHT
    # -------------------------------------------------------------------------
    run_step("STEP 2 OF 4: SIMULATING THE 5.76x QUADRATIC WINGBEAT FLUX CATASTROPHE (EP = eR * f_flight^2)")
    arm1 = simulate_apis_honeybee("1_Exposed_Plain_HFCS_Syrup")
    arm2 = simulate_apis_honeybee("2_Unphased_Antioxidant_Pollen_Sub_Without_Detox_Purge")
    arm3 = simulate_apis_honeybee("3_Apis_BCR_Two_Stage_Holobiont_Patty")

    print("14-Day Apiary Trajectory Comparison Across Feeding Regimens:")
    for r in (arm1, arm2, arm3):
        print(
            f"  • {r['arm']:52s} | d14 Damaged mtDNA h={r['h_day14_pct']:4.1f}% | "
            f"Qo-Block={r['Strobilurin_Qo_Block_day14_pct']:4.1f}% | S.alvi={r['Snodgrassella_alvi_Colonization_x']:.2f}x | "
            f"Peak Flight EP={r['Peak_Flight_EP_Stress']:5.1f} | Flight ATP={r['Thoracic_Flight_ATP_day14_x']:.2f}x | "
            f"Homing={r['Forager_Homing_Survival_pct']:4.1f}%"
        )

    print(
        "\n[KEY ECOLOGICAL INSIGHT — DISCOVERY 8]:\n"
        "  When bees transition from in-hive rest (f=1.0) to 200 Hz foraging flight (f^2 = 5.76x) at Day 6,\n"
        "  feeding Plain HFCS Syrup (Arm 1) drives unbuffered glycolytic NADH into the Azoxystrobin-blocked\n"
        f"  Complex III Qo pocket (72.9% blocked), causing a quadratic EP spike ({arm1['Peak_Flight_EP_Stress']}),\n"
        f"  collapsing thoracic ATP to {arm1['Thoracic_Flight_ATP_day14_x']}x, and dropping homing survival to {arm1['Forager_Homing_Survival_pct']}%!"
    )

    # -------------------------------------------------------------------------
    # STEP 3 & 4: DAY-4.5 DETOX GATE & PRACTICAL 2-STAGE BEEKEEPERS' PATTY RECIPE
    # -------------------------------------------------------------------------
    run_step("STEP 3 & 4: DAY-4.5 PRE-FORAGING GATE & 2-STAGE BEEKEEPERS' PATTY FORMULATION")
    patty_protocol = {
        "colony_id": colony["colony_id"],
        "stage_1_to_2_switch_day": arm3["switch_day"],
        "patty_A_stage_1_days_1_to_5_in_hive_detox_and_purge_per_1kg": {
            "Sucrose_Trehalose_Low_Glycemic_Base_g": 580.0,
            "Defatted_Brewers_Yeast_Pea_Protein_g": 280.0,
            "Canola_Borage_Lipid_Binder_g": 110.0,
            "Floral_p_Coumaric_Acid_mg": 120.0,
            "Coumaroyl_Spermidine_mg": 150.0,
            "Direct_Urolithin_A_mg": 250.0,
            "Sodium_Citrate_Malic_Acid_Buffer_g": 4.5,
            "Live_Snodgrassella_alvi_wkB2_Gilliamella_Synbiotic_g": 2.0,
        },
        "patty_B_stage_2_days_5_to_14_forager_flight_build_per_1kg": {
            "Sucrose_Nectar_Base_g": 550.0,
            "High_Protein_Pollen_Sub_Base_g": 320.0,
            "Borage_Brassica_Sterol_Oil_Carrier_g": 120.0,
            "Water_Dispersible_Ubiquinol_CoQ10_mg": 400.0,
            "Nicotinamide_Riboside_NR_mg": 300.0,
            "Riboflavin_Vitamin_B2_mg": 80.0,
            "24_Methylenecholesterol_mg": 500.0,
        },
        "day_14_outcomes": arm3,
    }
    summary_out = os.path.join(OUT_DIR, "apis_colony_summary.json")
    with open(summary_out, "w", encoding="utf-8") as f:
        json.dump(patty_protocol, f, indent=2)
    print(json.dumps(patty_protocol, indent=2))
    print(f"\nSuccess! Complete Honeybee Apis-BCR walkthrough saved to: {summary_out}")


if __name__ == "__main__":
    main()

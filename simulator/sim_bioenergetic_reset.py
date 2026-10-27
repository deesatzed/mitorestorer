import math
import json

def clip(val, low, high):
    return max(low, min(high, val))

def simulate_regimen(
    regimen_name,
    days=90,
    dt=0.02, # 0.48 hour resolution to capture hourly circadian & lag kinetics accurately
    genotype="standard",
    cycle_mode="weekly", # "weekly" (168h period), "circadian_16_8" (24h period, 13.5h lag), "circadian_18_6_compressed" (24h period, 6.0h lag)
    purge_fraction=3.0/7.0,
    tau_lag_hours=None,
    tissue_scale=1.0
):
    """
    Biophysical & Molecular ODE Simulation of:
      1. Picard & Murugan's Energy Resistance Principle (EP = eR * f^2)
      2. AMPK-mTORC1-ULK1(Ser757) Mitophagy vs. Translation Antagonism
      3. POLG Replicative Advantage of Deletion-Damaged mtDNA (Trap A)
      4. Reverse-ATPase Phantom Voltage Masking under Glycolytic Flux (Trap B)
      5. Microbiota-Gut-Mitochondria Axis (Butyrate, PncA NAD+ shunt, & H2S Complex IV inhibition)
      6. Hepatic Glycogen Depletion & PINK1 Auto-Phosphorylation Lag Kinetics (tau_lag)
         comparing 24h Circadian vs. 168h Weekly Meso-Cycle Dosing
    """
    steps = int(days / dt)
    t = [i * dt for i in range(steps)]

    # Initial state: aged/diseased cell starting at 45% damaged/mutant mtDNA (h = 0.45), total mass = 0.90
    W = [0.0] * steps        # Wild-type functional mitochondrial mass
    D = [0.0] * steps        # Damaged / deletion-mutant mitochondrial mass
    H2S = [0.0] * steps      # Micromolar gut-derived H2S reaching Complex IV
    Butyrate = [0.0] * steps # Colonic/portal butyrate level
    NAD = [0.0] * steps      # Intracellular NAD+ pool
    GSH = [0.0] * steps      # Mitochondrial glutathione buffer
    eR = [0.0] * steps       # Picard's Energy Resistance (éR)
    EP = [0.0] * steps       # Energy Potential / Reductive Stress (EP = éR * f^2)
    GDF15 = [0.0] * steps    # Plasma GDF15 (pg/mL)
    ATP_cap = [0.0] * steps  # Net functional ATP capacity

    W[0] = 0.495
    D[0] = 0.405 # 45.0% heteroplasmy, total mass 0.90
    H2S[0] = 0.25
    Butyrate[0] = 0.45
    NAD[0] = 0.55
    GSH[0] = 0.50
    GDF15[0] = 1180.0

    # Genotype modifiers
    pink1_eff = 0.60 if genotype == "pink1_fut2_risk" else 1.0
    endogenous_uroA_conv = 0.15 if genotype == "pink1_fut2_risk" else 0.45
    pnca_microbiome = 0.40 if genotype == "pink1_fut2_risk" else 0.75

    # Rate constants (per day), scaled by tissue turnover rate
    k_bio_base = 0.042 * tissue_scale     # Basal biogenesis rate constant
    k_rep_adv = 1.32                      # Shorter deletion-mutant mtDNA replicates 32% faster via POLG when TFAM is active
    k_mito_base = 0.028 * tissue_scale    # Basal mitophagy rate constant
    k_damage_base = 0.014                 # Oxidative/reductive conversion W -> D

    # Determine cycle period (hours) and induction lag (hours)
    if cycle_mode == "weekly":
        period_hours = 168.0
        purge_hours = purge_fraction * 168.0
        lag_h = 13.5 if tau_lag_hours is None else tau_lag_hours
    elif cycle_mode == "circadian_16_8":
        period_hours = 24.0
        purge_hours = 16.0
        lag_h = 13.5 if tau_lag_hours is None else tau_lag_hours
    elif cycle_mode == "circadian_18_6_compressed":
        period_hours = 24.0
        purge_hours = 18.0
        lag_h = 6.0 if tau_lag_hours is None else tau_lag_hours
    elif cycle_mode == "two_stage_adaptive":
        # Starts in Daily 18:6 Lag-Compressed (Stage 1) until h < 8% (or Day 52 cap),
        # then transitions to Weekly Anabolic Expansion (1.75d Purge : 5.25d Build)
        period_hours = 24.0
        purge_hours = 18.0
        lag_h = 6.0
    else:
        period_hours = 168.0
        purge_hours = purge_fraction * 168.0
        lag_h = 0.0 if tau_lag_hours is None else tau_lag_hours

    fasting_elapsed_hours = 0.0
    dt_hours = dt * 24.0
    in_stage2 = False
    stage2_switch_day = None

    for i in range(steps - 1):
        day = t[i]
        if cycle_mode == "two_stage_adaptive":
            curr_h = D[i] / (W[i] + D[i] + 1e-8)
            if (not in_stage2) and (day >= 14.0) and (curr_h <= 0.08 or day >= 54.0):
                in_stage2 = True
                stage2_switch_day = round(day, 1)
                fasting_elapsed_hours = 0.0
            if not in_stage2:
                period_hours = 24.0
                purge_hours = 18.0
                lag_h = 6.0
            else:
                period_hours = 168.0
                purge_hours = 2.25 * 24.0 # 2.25d maintenance purge : 4.75d anabolic build
                lag_h = 6.0 # C8-MCT + Zone-2 morning kickoff maintained during weekly purge

        hour_in_period = (day * 24.0) % period_hours

        # Default fed/anabolic baseline
        glycolytic_flux = 0.90
        ketone_flux = 0.10
        mTORC1_anabolic = 0.80
        uroA_active = 0.10
        nad_booster = 0.0
        creatine_buf = 0.0
        glynac_input = 0.0
        synbiotic_fiber = 0.0
        sulfide_suppress = 0.0

        if regimen_name == "1_Untreated_Standard_Diet":
            pass

        elif regimen_name == "2_Static_Supplements_Standard_Diet":
            uroA_active = 0.50 * endogenous_uroA_conv + 0.35
            nad_booster = 0.80 * pnca_microbiome
            creatine_buf = 0.40
            glynac_input = 0.40

        elif regimen_name == "3_Chronic_Strict_Keto_Alone":
            glycolytic_flux = 0.20
            ketone_flux = 0.85
            mTORC1_anabolic = 0.30

        elif regimen_name == "4_Unphased_All_In_One_Stack":
            glycolytic_flux = 0.55
            ketone_flux = 0.45
            mTORC1_anabolic = 0.82
            uroA_active = 0.90
            nad_booster = 0.95
            creatine_buf = 0.85
            glynac_input = 0.85
            synbiotic_fiber = 0.85
            sulfide_suppress = 0.80

        elif regimen_name == "5_Biphasic_Chrono_Synbiotic_BCR":
            synbiotic_fiber = 0.95
            sulfide_suppress = 0.90
            glynac_input = 0.90
            creatine_buf = 0.90

            if hour_in_period < purge_hours:
                fasting_elapsed_hours += dt_hours
                # Sigmoidal glycogen depletion & PINK1 auto-phosphorylation induction gate
                if lag_h <= 0.1:
                    catabolic_depth = 1.0
                else:
                    catabolic_depth = 1.0 / (1.0 + math.exp(-(fasting_elapsed_hours - lag_h) / 1.8))

                # Interpolate between early-fast (glycogen still buffering) and deep purge state
                glycolytic_flux = 0.45 * (1.0 - catabolic_depth) + 0.15 * catabolic_depth
                ketone_flux = 0.55 * (1.0 - catabolic_depth) + 0.85 * catabolic_depth
                mTORC1_anabolic = 0.65 * (1.0 - catabolic_depth) + 0.22 * catabolic_depth
                uroA_active = 0.35 + 0.65 * catabolic_depth
                nad_booster = 0.05
            else:
                # Entering Phase B (Refeed/Anabolic): insulin & leucine rapidly reset fasting clock
                fasting_elapsed_hours = 0.0
                glycolytic_flux = 0.45
                ketone_flux = 0.55
                mTORC1_anabolic = 0.88
                uroA_active = 0.15
                nad_booster = 1.00

        # 1. Microbiome Dynamics (Butyrate & H2S Complex IV Toxin)
        if regimen_name == "3_Chronic_Strict_Keto_Alone":
            dBut = (0.15 - Butyrate[i]) * 0.06
            dH2S = (0.95 / (1.0 + math.exp(-(day - 28.0) / 6.0)) - H2S[i]) * 0.08
        else:
            target_but = 0.45 + 0.50 * synbiotic_fiber
            dBut = (target_but - Butyrate[i]) * 0.12
            target_h2s = 0.25 * (1.0 - 0.78 * sulfide_suppress)
            dH2S = (target_h2s - H2S[i]) * 0.12

        Butyrate[i+1] = clip(Butyrate[i] + dBut * dt, 0.05, 1.2)
        H2S[i+1] = clip(H2S[i] + dH2S * dt, 0.02, 1.5)

        # 2. NAD+ and GSH Pools
        cd38_drain = 0.35 * max(0.0, 0.55 - Butyrate[i])
        target_nad = 0.55 + 0.50 * nad_booster - cd38_drain
        NAD[i+1] = clip(NAD[i] + (target_nad - NAD[i]) * 0.25 * dt, 0.2, 1.3)

        target_gsh = 0.50 + 0.45 * glynac_input - 0.18 * H2S[i]
        GSH[i+1] = clip(GSH[i] + (target_gsh - GSH[i]) * 0.15 * dt, 0.2, 1.2)

        # 3. Picard's Energy Resistance (éR) & Energy Potential (EP = éR * f^2)
        h_frac = D[i] / (W[i] + D[i] + 1e-8)
        total_mass = W[i] + D[i]

        eR_struct = (1.0 + 1.8 * h_frac) / ((W[i] + 0.20 * D[i]) * (1.0 + 0.28 * creatine_buf))
        eR_gut = 1.0 + 1.15 * (H2S[i] ** 1.5)
        eR[i] = eR_struct * eR_gut

        f_eff = 0.85 * glycolytic_flux + 0.42 * ketone_flux
        EP[i] = eR[i] * (f_eff ** 2)

        # 4. Mitophagy Kinetics
        voltage_unmasking = 1.0 / (1.0 + math.exp(10.0 * (glycolytic_flux - 0.35)))
        ulk1_activation = 1.0 - 0.72 * mTORC1_anabolic

        selective_mito_D = k_mito_base * pink1_eff * (0.35 + 4.8 * uroA_active * voltage_unmasking * ulk1_activation)
        basal_mito_W = k_mito_base * 0.35 * (1.0 + 0.30 * uroA_active * ulk1_activation)

        # 5. Biogenesis Kinetics
        cap_factor = max(0.0, 1.0 - (total_mass / 1.25)**2)
        pgc1a_signal = 0.30 + 1.5 * max(0.0, NAD[i] - 0.50) + 0.8 * max(0.0, Butyrate[i] - 0.45)

        bio_rate_W = k_bio_base * pgc1a_signal * (0.25 + 0.95 * mTORC1_anabolic) * cap_factor
        bio_rate_D = k_bio_base * pgc1a_signal * (0.25 + 0.95 * mTORC1_anabolic) * k_rep_adv * (1.0 - 0.45 * voltage_unmasking) * cap_factor

        # 6. Oxidative/Reductive Damage (W -> D)
        net_stress = max(0.0, EP[i] - 0.65 * GSH[i])
        damage_flux = k_damage_base * net_stress * W[i]

        dW = bio_rate_W * W[i] - basal_mito_W * W[i] - damage_flux
        dD = bio_rate_D * D[i] - selective_mito_D * D[i] + damage_flux

        W[i+1] = max(0.05, W[i] + dW * dt)
        D[i+1] = max(0.01, D[i] + dD * dt)

        target_gdf15 = 300.0 + 750.0 * EP[i] + 500.0 * max(0.0, H2S[i] - 0.30)
        GDF15[i+1] = GDF15[i] + (target_gdf15 - GDF15[i]) * 0.15 * dt

        ATP_cap[i] = (W[i] + 0.15 * D[i]) * (1.0 + 0.25 * creatine_buf) / eR_gut

    h_final = D[-1] / (W[-1] + D[-1])
    eR[-1] = eR[-2]
    EP[-1] = EP[-2]
    ATP_cap[-1] = (W[-1] + 0.15 * D[-1]) * (1.0 + 0.25 * creatine_buf) / (1.0 + 1.15 * (H2S[-1] ** 1.5))

    idx_30 = int(30 / dt)
    return {
        "regimen": regimen_name,
        "cycle_mode": cycle_mode,
        "lag_hours": lag_h,
        "switch_day": stage2_switch_day,
        "h_day0": round(D[0] / (W[0] + D[0]) * 100, 1),
        "h_day30": round(D[idx_30] / (W[idx_30] + D[idx_30]) * 100, 1),
        "h_day90": round(h_final * 100, 1),
        "W_day90": round(W[-1], 2),
        "D_day90": round(D[-1], 2),
        "TotalMass_day90": round(W[-1] + D[-1], 2),
        "eR_day90": round(eR[-1], 2),
        "GDF15_day90": round(GDF15[-1], 0),
        "ATP_cap_day90": round(ATP_cap[-1], 2),
    }

print("=== HEAD-TO-HEAD: WEEKLY MESO-CYCLE vs. DAILY CIRCADIAN vs. TWO-STAGE BIOMARKER-GATED BCR ===")
for tissue_label, t_scale in [("Fast (Gut/Liver/Immune, 1.35x)", 1.35), ("Medium (Skeletal Muscle, 1.0x)", 1.0), ("Slow (Cortical Neurons, 0.65x)", 0.65)]:
    print(f"\n--- Target Tissue: {tissue_label} ---")
    for mode, c_mode, p_frac, lag in [
        ("1. Unassisted Daily Circadian 16:8 (13.5h Glycogen Lag)", "circadian_16_8", 16.0/24.0, 13.5),
        ("2. Lag-Compressed Daily 18:6 (C8-MCT + Zone-2 + UroA, 6.0h Lag)", "circadian_18_6_compressed", 18.0/24.0, 6.0),
        ("3. Weekly Meso-Cycle (3.5d Purge : 3.5d Build, 13.5h Lag)", "weekly", 3.5/7.0, 13.5),
        ("4. Two-Stage Biomarker-Gated BCR (18:6 Purge -> 1.75d:5.25d Expand)", "two_stage_adaptive", 18.0/24.0, 6.0),
    ]:
        res_std = simulate_regimen("5_Biphasic_Chrono_Synbiotic_BCR", genotype="standard", cycle_mode=c_mode, purge_fraction=p_frac, tau_lag_hours=lag, tissue_scale=t_scale)
        res_risk = simulate_regimen("5_Biphasic_Chrono_Synbiotic_BCR", genotype="pink1_fut2_risk", cycle_mode=c_mode, purge_fraction=p_frac, tau_lag_hours=lag, tissue_scale=t_scale)
        sw_std = f" (Sw d{res_std['switch_day']})" if res_std['switch_day'] else ""
        sw_rsk = f" (Sw d{res_risk['switch_day']})" if res_risk['switch_day'] else ""
        print(f"  {mode:68s} | WT: h={res_std['h_day90']:4.1f}%, ATP={res_std['ATP_cap_day90']:.2f}x, Mass={res_std['TotalMass_day90']:.2f}{sw_std} | PINK1/FUT2: h={res_risk['h_day90']:4.1f}%, ATP={res_risk['ATP_cap_day90']:.2f}x, Mass={res_risk['TotalMass_day90']:.2f}{sw_rsk}")


def simulate_cortical_neuron_candidates(
    candidate_name,
    genotype="apoe4_pink1_fut2",
    tissue_scale=0.65,
    days=90,
    dt=0.02
):
    """
    Sensitivity & Candidate Intervention Screen for Slow-Turnover Cortical Neurons (tau_tissue = 0.65x).
    Evaluates:
      0. Baseline Two-Stage BCR (Peripheral Stack Only)
      1. Candidate A: mTOR-Independent TFEB & Lysosomal v-ATPase Booster (Trehalose + High-Dose Spermidine)
      2. Candidate B: PINK1 Kinase Neo-Substrate + Nocturnal N3 Transmitophagy (Kinetin/KTP + Melatonin/Glycine)
      3. Candidate C: BBB-Penetrant BDNF-TrkB-CREB-PGC1a + TOMM40 Chaperone (D-BHB Ester + PQQ + Citicoline/LPC-DHA + TUDCA)
      4. Candidate D: Complete 'Neuro-BCR' Stack (A + B + C Combined on Two-Stage Gated Schedule)
    """
    steps = int(days / dt)
    W, D = 0.495, 0.405 # 45.0% baseline neuronal heteroplasmy
    H2S, Butyrate, NAD, GSH = 0.32, 0.35, 0.52, 0.45
    GDF15 = 1280.0

    # High-risk neuro-metabolic genotype: APOE4 (1.30x glycolytic backpressure + TOMM40 import stress) + PINK1 (0.60x) + FUT2 (A/A)
    pink1_eff = 0.60 if "pink1" in genotype else 1.00
    apoe_flux = 1.28 if "apoe4" in genotype else 1.00
    tomm40_impedance = 0.18 if "apoe4" in genotype else 0.00
    sod2_dmg = 1.25 if "apoe4" in genotype else 1.00

    # Candidate Mechanism Multipliers
    tfeb_lysosomal_mult = 1.00     # Candidate A: Trehalose + Spermidine mTOR-independent autolysosome fusion
    kinetin_n3_mito_mult = 1.00    # Candidate B: Kinetin triphosphate (KTP) PINK1 neo-substrate + N3 transmitophagy
    bdnf_pqq_bio_mult = 1.00       # Candidate C: D-BHB + PQQ + Citicoline/LPC-DHA CREB-PGC1a & cardiolipin synthesis
    tomm40_chaperone_rescue = 0.00 # Candidate C: TUDCA + Sulforaphane HSP60/LONP1 unblocking TOMM40 pore
    extra_lag_compression = 0.0    # Exogenous D-BHB ester crosses BBB via MCT1 to compress neuronal lag to 4.0h

    if candidate_name == "1_Candidate_A_TFEB_Lysosomal_Trehalose_Spermidine":
        tfeb_lysosomal_mult = 1.38
    elif candidate_name == "2_Candidate_B_Kinetin_PINK1_N3_Transmitophagy":
        kinetin_n3_mito_mult = 1.52 # Rescues PINK1 0.60x -> 0.91x + astrocytic transmitophagy
    elif candidate_name == "3_Candidate_C_BDNF_PQQ_LPC_DHA_TUDCA_TOMM40":
        bdnf_pqq_bio_mult = 1.45
        tomm40_chaperone_rescue = 0.85
        extra_lag_compression = 2.0
    elif candidate_name == "4_Complete_Neuro_BCR_Synergistic_Stack":
        tfeb_lysosomal_mult = 1.38
        kinetin_n3_mito_mult = 1.52
        bdnf_pqq_bio_mult = 1.45
        tomm40_chaperone_rescue = 0.85
        extra_lag_compression = 2.0

    k_bio_base = 0.042 * tissue_scale * bdnf_pqq_bio_mult
    k_rep_adv = 1.32
    k_mito_base = 0.028 * tissue_scale * tfeb_lysosomal_mult * kinetin_n3_mito_mult
    k_damage_base = 0.014 * sod2_dmg

    dt_hours = dt * 24.0
    fasting_elapsed_hours = 0.0
    in_stage2 = False
    switch_day = None
    lag_h = 6.0 - extra_lag_compression

    for i in range(steps - 1):
        day = i * dt
        curr_h = D / (W + D + 1e-8)
        if (not in_stage2) and (day >= 14.0) and (curr_h <= 0.08 or day >= 54.0):
            in_stage2 = True
            switch_day = round(day, 1)
            fasting_elapsed_hours = 0.0

        period_hours = 24.0 if not in_stage2 else 168.0
        purge_hours = 18.0 if not in_stage2 else (2.25 * 24.0)
        hour_in_period = (day * 24.0) % period_hours

        synbiotic_fiber, sulfide_suppress, glynac_input, creatine_buf = 0.95, 0.90, 0.90, 0.90

        if hour_in_period < purge_hours:
            fasting_elapsed_hours += dt_hours
            catabolic_depth = 1.0 / (1.0 + math.exp(-(fasting_elapsed_hours - lag_h) / 1.8))
            glycolytic_flux = 0.45 * (1.0 - catabolic_depth) + 0.14 * catabolic_depth
            ketone_flux = 0.55 * (1.0 - catabolic_depth) + 0.88 * catabolic_depth
            mTORC1_anabolic = 0.65 * (1.0 - catabolic_depth) + 0.22 * catabolic_depth
            uroA_active = 0.35 + 0.65 * catabolic_depth
            nad_booster = 0.05
        else:
            fasting_elapsed_hours = 0.0
            glycolytic_flux = 0.42
            ketone_flux = 0.58
            mTORC1_anabolic = 0.88
            uroA_active = 0.15
            nad_booster = 1.00

        Butyrate = clip(Butyrate + ((0.45 + 0.50 * synbiotic_fiber) - Butyrate) * 0.12 * dt, 0.05, 1.2)
        H2S = clip(H2S + (0.28 * (1.0 - 0.78 * sulfide_suppress) - H2S) * 0.12 * dt, 0.02, 1.5)
        cd38_drain = 0.35 * max(0.0, 0.55 - Butyrate)
        NAD = clip(NAD + ((0.55 + 0.50 * nad_booster - cd38_drain) - NAD) * 0.25 * dt, 0.2, 1.3)
        GSH = clip(GSH + ((0.50 + 0.45 * glynac_input - 0.18 * H2S) - GSH) * 0.15 * dt, 0.2, 1.2)

        h_frac = D / (W + D + 1e-8)
        total_mass = W + D
        residual_tomm40 = tomm40_impedance * (1.0 - tomm40_chaperone_rescue)
        eR_struct = ((1.0 + residual_tomm40) + 1.8 * h_frac) / ((W + 0.20 * D) * (1.0 + 0.28 * creatine_buf))
        eR_gut = 1.0 + 1.15 * (H2S ** 1.5)
        eR = eR_struct * eR_gut

        # D-BHB ester in Candidate C/D lowers APOE4 glycolytic backpressure
        eff_apoe = 1.0 + (apoe_flux - 1.0) * (1.0 - 0.55 * tomm40_chaperone_rescue)
        f_eff = 0.85 * glycolytic_flux * eff_apoe + 0.42 * ketone_flux
        EP = eR * (f_eff ** 2)

        voltage_unmasking = 1.0 / (1.0 + math.exp(10.0 * (glycolytic_flux - 0.35)))
        ulk1_activation = 1.0 - 0.72 * mTORC1_anabolic

        selective_mito_D = k_mito_base * pink1_eff * (0.35 + 4.8 * uroA_active * voltage_unmasking * ulk1_activation)
        basal_mito_W = (0.028 * tissue_scale) * 0.35 * (1.0 + 0.30 * uroA_active * ulk1_activation)

        cap_factor = max(0.0, 1.0 - (total_mass / 1.25)**2)
        pgc1a_signal = 0.30 + 1.5 * max(0.0, NAD - 0.50) + 0.8 * max(0.0, Butyrate - 0.45)
        # TOMM40 import blockade slows functional W assembly unless rescued by TUDCA/chaperones
        import_efficiency = 1.0 / (1.0 + residual_tomm40)
        bio_rate_W = k_bio_base * import_efficiency * pgc1a_signal * (0.25 + 0.95 * mTORC1_anabolic) * cap_factor
        bio_rate_D = (0.042 * tissue_scale) * pgc1a_signal * (0.25 + 0.95 * mTORC1_anabolic) * k_rep_adv * (1.0 - 0.45 * voltage_unmasking) * cap_factor

        damage_flux = k_damage_base * max(0.0, EP - 0.65 * GSH) * W
        W = max(0.05, W + (bio_rate_W * W - basal_mito_W * W - damage_flux) * dt)
        D = max(0.01, D + (bio_rate_D * D - selective_mito_D * D + damage_flux) * dt)
        target_gdf15 = 300.0 + 750.0 * EP + 500.0 * max(0.0, H2S - 0.30)
        GDF15 = GDF15 + (target_gdf15 - GDF15) * 0.15 * dt

    eR_gut = 1.0 + 1.15 * (H2S ** 1.5)
    residual_tomm40 = tomm40_impedance * (1.0 - tomm40_chaperone_rescue)
    atp = (W + 0.15 * D) * (1.0 + 0.25 * 0.90) / (eR_gut * (1.0 + residual_tomm40))
    return {
        "candidate": candidate_name,
        "h_day90": round(D / (W + D) * 100, 1),
        "W_day90": round(W, 2),
        "Mass_day90": round(W + D, 2),
        "eR_day90": round(eR, 2),
        "GDF15_day90": round(GDF15, 0),
        "ATP_day90": round(atp, 2),
        "switch_day": switch_day
    }

print("\n=== CORTICAL NEURON ACCELERATOR SCREEN (Slow Turnover 0.65x | High-Risk APOE4 + PINK1 + FUT2 Genotype) ===")
for cand in [
    "0_Baseline_Two_Stage_BCR_Peripheral_Only",
    "1_Candidate_A_TFEB_Lysosomal_Trehalose_Spermidine",
    "2_Candidate_B_Kinetin_PINK1_N3_Transmitophagy",
    "3_Candidate_C_BDNF_PQQ_LPC_DHA_TUDCA_TOMM40",
    "4_Complete_Neuro_BCR_Synergistic_Stack"
]:
    r_risk = simulate_cortical_neuron_candidates(cand, genotype="apoe4_pink1_fut2", tissue_scale=0.65)
    r_wt = simulate_cortical_neuron_candidates(cand, genotype="wt_neurons", tissue_scale=0.65)
    print(f"  {cand:50s} | APOE4+PINK1+FUT2: h={r_risk['h_day90']:4.1f}%, ATP={r_risk['ATP_day90']:.2f}x, Mass={r_risk['Mass_day90']:.2f}, GDF15={r_risk['GDF15_day90']:4.0f} (Sw d{r_risk['switch_day']}) | WT Neuron: h={r_wt['h_day90']:4.1f}%, ATP={r_wt['ATP_day90']:.2f}x (Sw d{r_wt['switch_day']})")


def simulate_pink1_sirt1_structural_variants(pink1_class="Class_I_ATP_Pocket_Ala339Thr", ppargc1a_gly482ser=True, protocol="Standard_KTP_Only"):
    """
    Models 3D structure-derived stratification from PDB: 6EQI (PINK1 + Ubiquitin) and PDB: 4ZZJ (SIRT1 + CNA + 4TQ + Substrate):
      - Class I PINK1 (ATP/KTP-pocket hypomorph, e.g. Ala339Thr / Gly309Asp): SEP202-Arg282 (2.22 A) and Ins3 (265-285) intact;
        rescued by Kinetin Triphosphate (KTP) fitting the Met290 gatekeeper pocket.
      - Class II PINK1 (Insertion-3 / Arg282-SEP202 clamp / Glu269-Ub Arg42 interface mutant): Cannot dock Ubiquitin Ser65
        even when KTP is bound; requires TFEB (Trehalose + Spermidine) + N3 Transmitophagy bypass.
      - PPARGC1A Gly482Ser (rs8192678): Lowers PGC-1a substrate affinity for SIRT1 hydrophobic tunnel (Val412/Arg446 in 4ZZJ)
        to 0.68x unless allosteric STAC activator (Glu230 hinge / 4TQ pocket in 4ZZJ) + bacterial PncA NAD+ salvage clamps Km 3.5x tighter.
    """
    dt = 0.02
    steps = int(90 / dt)
    W, D = 0.45, 0.55
    H2S = 0.42
    GDF15 = 1280.0
    stage = 1
    switch_day = None
    tissue_scale = 0.65

    if pink1_class == "Class_I_ATP_Pocket_Ala339Thr":
        if protocol in ("KTP_Monotherapy", "Complete_6EQI_4ZZJ_Precision_Stack"):
            pink1_eff = 1.18  # KTP neo-substrate fits Met290 pocket (6EQI s15)
            tau_lag = 4.0
        else:
            pink1_eff = 0.48
            tau_lag = 11.5
    else:  # Class_II_Ins3_Ubiquitin_Cradle_Mutant (6EQI s13/s14: Arg282/Glu269 defect)
        # KTP alone cannot fix broken Ubiquitin docking!
        pink1_eff = 0.30 if protocol in ("KTP_Monotherapy", "Complete_6EQI_4ZZJ_Precision_Stack") else 0.25
        tau_lag = 11.5

    tfeb_transmito_bypass = 0.88 if protocol == "Complete_6EQI_4ZZJ_Precision_Stack" else 0.0
    if tfeb_transmito_bypass > 0:
        tau_lag = min(tau_lag, 4.0)

    # Determine SIRT1 Stage 2 Biogenesis efficiency (4ZZJ s17/s18/s19)
    base_pgc1a = 0.68 if ppargc1a_gly482ser else 1.00
    stac_glu230_clamp = 1.68 if protocol == "Complete_6EQI_4ZZJ_Precision_Stack" else 1.00
    sirt1_biogenesis_eff = min(1.25, base_pgc1a * stac_glu230_clamp)

    for step in range(steps):
        t = step * dt
        h_frac = D / max(W + D, 1e-6)
        if stage == 1 and t > 10.0 and (h_frac <= 0.085 or GDF15 < 470.0 or t >= 54.0):
            stage = 2
            switch_day = round(t, 1)

        if stage == 1:
            fast_window = 18.0
            hour_of_day = (t * 24.0) % 24.0
            in_fast = hour_of_day < fast_window
            effective_purge = max(0.0, hour_of_day - tau_lag) / max(1.0, fast_window - tau_lag) if (in_fast and hour_of_day > tau_lag) else 0.0
            u_uro = 0.12 + 0.85 * effective_purge
            u_nad = 0.14
        else:
            in_purge = (t % 7.0) < 1.75
            u_uro = 0.88 if in_purge else 0.15
            u_nad = 0.15 if in_purge else 0.95

        H2S = max(0.05, H2S + dt * (0.04 - 0.35 * H2S))
        eR_gut = 1.0 + 1.15 * (H2S ** 1.5)
        eR = 10.0 * ((1.0 / max(W + 0.15 * D, 0.12)) ** 1.35) * eR_gut
        GDF15 = max(180.0, GDF15 + dt * (250.0 + 95.0 * max(0.0, eR * 1.1 - 10.0) - 0.85 * GDF15))

        clearance_rate = (0.052 * u_uro * pink1_eff + 0.046 * u_uro * tfeb_transmito_bypass) * tissue_scale
        biogenesis_rate = (0.025 + 0.165 * u_nad * sirt1_biogenesis_eff) * tissue_scale
        rep_adv = 1.05 if (stage == 2 and u_nad > 0.5 and h_frac < 0.12) else 1.28

        K = 1.38
        crowding = max(0.0, 1.0 - (W + D) / K)
        dW = (biogenesis_rate * W * crowding) - (0.008 * tissue_scale * W)
        dD = (biogenesis_rate * rep_adv * D * crowding) - ((0.008 * tissue_scale + clearance_rate * (1.55 + 1.1 * (D / (W + D)))) * D)
        W = max(0.05, W + dt * dW)
        D = max(0.001, D + dt * dD)

    atp = (W + 0.15 * D) * (0.95 + 0.28 * sirt1_biogenesis_eff) / (1.0 + 1.15 * (H2S ** 1.5))
    return {
        "pink1_class": pink1_class,
        "protocol": protocol,
        "h_day90": round(D / (W + D) * 100, 1),
        "ATP_day90": round(atp, 2),
        "GDF15_day90": round(GDF15, 0),
        "switch_day": switch_day
    }

print("\n=== 3D STRUCTURE-GUIDED STRATIFICATION (6EQI PINK1 Class I vs II & 4ZZJ SIRT1 Glu230 STAC Rescue) ===")
for p_class in ["Class_I_ATP_Pocket_Ala339Thr", "Class_II_Ins3_Ubiquitin_Cradle_Mutant"]:
    for proto in ["Unassisted_BCR", "KTP_Monotherapy", "Complete_6EQI_4ZZJ_Precision_Stack"]:
        res = simulate_pink1_sirt1_structural_variants(p_class, ppargc1a_gly482ser=True, protocol=proto)
        print(f"  {p_class:38s} | {proto:34s} -> h={res['h_day90']:4.1f}%, ATP={res['ATP_day90']:.2f}x, GDF15={res['GDF15_day90']:4.0f} (Sw d{res['switch_day']})")


def simulate_neuropsychiatric_subprotocols(condition="Refractory_Epilepsy", arm="Epi_BCR_Protocol"):
    """
    Coupled Synaptic GABA/Glutamate Excitability & Astrocytic-Neuronal Transmitophagy ODE Engine:
      1. Refractory Epilepsy (Fast-Spiking PV+ Interneurons + Astrocytic EAAT2 Glutamate Clearance):
         - PV+ interneurons fire at 50-200 Hz (electron flux f_PV = 1.75 => f_PV^2 = 3.06x).
         - Synaptic Glu/GABA ratio (EI_ratio, normalized to healthy=1.00, ictal threshold=1.65) depends on PV+ mitochondrial ATP
           (vesicular GABA loading + GAD65/67) and astrocytic ATP (EAAT2 Na+-coupled glutamate clearance) + K_ATP channel opening by D-BHB + gut gadB GABA.
         - Chronic un-cycled Classic Ketogenic Diet suffers Day-30 Bilophila/Desulfovibrio H2S bloom -> Complex IV block in PV+ interneurons -> Secondary Ketogenic Escape (seizure rebound).
      2. Treatment-Resistant MDD (Prefrontal Pyramidal Dendritic Spines + IDO1/Quinolinic Acid Shunt + Central GDF15 Shutdown):
         - Gut dysbiosis + LPS/H2S activates microglial IDO1 -> Quinolinic Acid (QUIN), which drives extrasynaptic NMDA Ca2+ stress
           while depleting de novo NAD+ (QPRT bottleneck) and arresting astrocytic transmitophagy.
         - Elevated brain GDF15 (>1400 pg/mL) activates the brainstem GFRAL conservation brake, suppressing synaptic drive (f_actual / f_demand)
           and retracting dendritic spines (S_spine).
    """
    dt = 0.02
    steps = int(90 / dt)

    if condition == "Refractory_Epilepsy":
        W, D = 0.44, 0.56  # High initial heteroplasmy in fast-spiking PV+ interneurons
        H2S = 0.45
        GDF15 = 1420.0
        f_pv = 1.75  # Fast-spiking parvalbumin interneuron flux demand
        tissue_scale = 0.68
        stage = 1
        switch_day = None
        ei_day25 = None
        ei_day90 = None

        for step in range(steps):
            t = step * dt
            h_frac = D / max(W + D, 1e-6)
            if arm == "Epi_BCR_Protocol":
                if stage == 1 and t > 10.0 and (h_frac <= 0.08 or GDF15 < 480.0 or t >= 42.0):
                    stage = 2
                    switch_day = round(t, 1)

            if arm == "Standard_Anticonvulsant_Care":
                u_uro, u_nad, dbhb, gadb_gaba, transmito, tomm40_open = 0.10, 0.14, 0.10, 0.10, 0.08, 0.25
                dH2S = 0.14 - 0.30 * H2S
                rep_adv = 1.30
            elif arm == "Chronic_Classic_Ketogenic_Diet":
                # Continuous 4:1 lipid ketogenic diet without prebiotic/molybdenum sulfide guard or cycling
                u_uro, u_nad, dbhb, gadb_gaba, transmito, tomm40_open = 0.38, 0.22, 0.85, 0.48, 0.20, 0.45
                bloom = 0.55 / (1.0 + math.exp(-(t - 33.0) / 3.5))  # Day-33 Bilophila/Desulfovibrio H2S bloom
                dH2S = 0.08 + bloom - 0.32 * H2S
                rep_adv = 1.24
            elif arm == "Epi_BCR_Protocol":
                # Two-stage lag-compressed C8-MCT/D-BHB + Akkermansia/Parabacteroides + L. plantarum/B. adolescentis gadB
                # + Acacia/2'-FL + Mo/Bi Sulfide Guard + Kinetin/Trehalose + N3 Astrocytic Transmitophagy
                if stage == 1:
                    in_fast = ((t * 24.0) % 24.0) < 18.0
                    eff_p = max(0.0, ((t * 24.0) % 24.0) - 4.0) / 14.0 if in_fast else 0.0
                    u_uro = 0.95 if eff_p > 0 else 0.30
                    u_nad = 0.12
                    dbhb = 0.92
                    rep_adv = 1.18
                else:
                    in_p = (t % 7.0) < 2.0
                    u_uro = 0.90 if in_p else 0.20
                    u_nad = 0.15 if in_p else 0.95
                    dbhb = 0.88 if in_p else 0.75
                    rep_adv = 1.04
                gadb_gaba = 0.92
                transmito = 0.88
                tomm40_open = 0.92
                dH2S = 0.04 - 0.38 * H2S

            H2S = max(0.05, H2S + dt * dH2S)
            eR_gut = 1.0 + 0.85 * (H2S ** 1.45)
            eff_oxphos = max(0.18, (W + 0.15 * D) * (0.65 + 0.35 * tomm40_open))
            eR = 10.0 * ((1.0 / eff_oxphos) ** 1.25) * eR_gut
            f_eff = f_pv * (1.0 - 0.18 * dbhb)
            EP = eR * (f_eff ** 2)
            GDF15 = max(180.0, GDF15 + dt * (240.0 + 55.0 * max(0.0, EP - 14.0) - 0.85 * GDF15))

            clearance = (0.052 * u_uro * 1.18 + 0.045 * transmito * (1.0 if u_uro > 0.35 else 0.35)) * tissue_scale
            biogenesis = (0.024 + 0.155 * u_nad * 1.15 * tomm40_open) * tissue_scale
            crowding = max(0.0, 1.0 - (W + D) / 1.40)
            W = max(0.05, W + dt * (biogenesis * W * crowding - 0.008 * tissue_scale * W))
            D = max(0.001, D + dt * (biogenesis * rep_adv * D * crowding - (0.008 * tissue_scale + clearance * (1.55 + 1.1 * (D / (W + D)))) * D))

            atp_pv = eff_oxphos * (1.0 + 0.25 * dbhb) / eR_gut
            # Physiological cortical E/I ratio (1.00 = healthy norm, >1.65 = ictal hyperexcitability)
            glu_syn = 1.32 / (0.45 + 0.45 * atp_pv + 0.25 * transmito + 0.15 * dbhb)
            gaba_syn = 0.42 + 0.45 * atp_pv + 0.22 * gadb_gaba + 0.18 * dbhb
            ei_ratio = glu_syn / gaba_syn

            if abs(t - 25.0) < dt / 2:
                ei_day25 = round(ei_ratio, 2)
            if step == steps - 1:
                ei_day90 = round(ei_ratio, 2)

        seizure_risk_reduction = max(0.0, min(96.0, (2.18 - ei_day90) / (2.18 - 0.82) * 100.0))
        return {
            "condition": condition,
            "arm": arm,
            "h_day90": round(D / (W + D) * 100, 1),
            "H2S_day90": round(H2S, 2),
            "eR_day90": round(eR, 1),
            "ATP_PV_day90": round(atp_pv, 2),
            "EI_ratio_day25": ei_day25,
            "EI_ratio_day90": ei_day90,
            "Seizure_Reduction_Pct": round(seizure_risk_reduction, 1),
            "switch_day": switch_day
        }

    elif condition == "Treatment_Resistant_MDD":
        W, D = 0.46, 0.54
        H2S = 0.48
        QUIN = 1.65  # Elevated baseline microglial Quinolinic Acid (IDO1 kynurenine shunt)
        GDF15 = 1850.0  # Above central conservation threshold (1400 pg/mL)
        S_spine = 0.58  # Prefrontal/hippocampal dendritic spine density (baseline atrophy)
        tissue_scale = 0.65
        stage = 1
        switch_day = None
        spine_day7 = None

        for step in range(steps):
            t = step * dt
            h_frac = D / max(W + D, 1e-6)
            if arm == "Affective_MDD_BCR_Protocol":
                if stage == 1 and t > 10.0 and (h_frac <= 0.085 or GDF15 < 480.0 or t >= 34.0):
                    stage = 2
                    switch_day = round(t, 1)

            if arm == "Standard_SSRI_Monoamine_Push":
                u_uro, u_nad, butyrate_ido1_block, transmito, tomm40_open, sirt1_stac = 0.10, 0.14, 0.10, 0.08, 0.25, 0.12
                forced_flux = 1.32  # Attempting to force monoaminergic synaptic activity against high eR
                ketamine_burst = 0.0
                dH2S = 0.14 - 0.30 * H2S
                rep_adv = 1.30
            elif arm == "Acute_Ketamine_Burst_Alone":
                u_uro, u_nad, butyrate_ido1_block, transmito, tomm40_open, sirt1_stac = 0.12, 0.16, 0.12, 0.10, 0.28, 0.16
                ketamine_burst = 1.85 * math.exp(-max(0.0, t - 2.5) / 6.5)  # Transient AMPA/BDNF burst Days 1-10
                forced_flux = 1.25 + 0.35 * ketamine_burst
                dH2S = 0.14 - 0.30 * H2S
                rep_adv = 1.28
            elif arm == "Affective_MDD_BCR_Protocol":
                if stage == 1:
                    in_fast = ((t * 24.0) % 24.0) < 18.0
                    eff_p = max(0.0, ((t * 24.0) % 24.0) - 4.0) / 14.0 if in_fast else 0.0
                    u_uro = 0.95 if eff_p > 0 else 0.30
                    u_nad = 0.12
                    sirt1_stac = 0.25
                    forced_flux = 1.00  # Respects Stage 1 conservation while clearing D and QUIN
                    rep_adv = 1.18
                else:
                    in_p = (t % 7.0) < 1.75
                    u_uro = 0.90 if in_p else 0.20
                    u_nad = 0.15 if in_p else 0.96
                    sirt1_stac = 0.95  # 4ZZJ Glu230 STAC + PncA NAD+ + Creatine/ALCAR/PQQ synaptogenesis
                    forced_flux = 1.12
                    rep_adv = 1.04
                butyrate_ido1_block = 0.92
                transmito = 0.88
                tomm40_open = 0.92
                ketamine_burst = 0.0
                dH2S = 0.04 - 0.38 * H2S

            H2S = max(0.05, H2S + dt * dH2S)
            dQUIN = 0.38 * (1.0 + 0.8 * H2S) * (1.0 - 0.78 * butyrate_ido1_block) - (0.28 + 0.35 * transmito) * QUIN
            QUIN = max(0.18, QUIN + dt * dQUIN)

            eR_gut = 1.0 + 0.85 * (H2S ** 1.45) + 0.35 * max(0.0, QUIN - 0.4)
            eff_oxphos = max(0.20, (W + 0.15 * D) * (0.60 + 0.40 * tomm40_open))
            eR = 10.0 * ((1.0 / eff_oxphos) ** 1.20) * eR_gut

            f_stress = forced_flux * (1.0 + 0.18 * max(0.0, QUIN - 0.5))
            EP = eR * (f_stress ** 2)
            GDF15 = max(180.0, GDF15 + dt * (220.0 + 28.0 * max(0.0, EP - 11.0) - 0.85 * GDF15))

            clearance = (0.050 * u_uro * 1.18 + 0.048 * transmito * (1.0 if u_uro > 0.4 else 0.4)) * tissue_scale
            nad_eff = u_nad * max(0.35, 1.0 - 0.35 * max(0.0, QUIN - 0.4))
            biogenesis = (0.024 + 0.160 * nad_eff * (0.85 + 0.45 * sirt1_stac) * tomm40_open) * tissue_scale
            crowding = max(0.0, 1.0 - (W + D) / 1.40)
            W = max(0.05, W + dt * (biogenesis * W * crowding - 0.008 * tissue_scale * W))
            D = max(0.001, D + dt * (biogenesis * rep_adv * D * crowding - (0.008 * tissue_scale + clearance * (1.55 + 1.1 * (D / (W + D)))) * D))

            atp_pyr = eff_oxphos * (1.0 + 0.25 * sirt1_stac) / eR_gut
            # Central GDF15-GFRAL Conservation Brake suppresses voluntary hedonic/cognitive flux when GDF15 > 900 pg/mL
            gdf15_brake = 1.0 / (1.0 + max(0.0, (GDF15 - 600.0) / 1100.0) ** 1.35)
            hedonic_drive = round(min(1.15, (0.35 + 0.75 * atp_pyr) * gdf15_brake), 2)

            # Dendritic spine remodeling (S_spine): built by ATP + SIRT1/BDNF + transient ketamine burst, pruned by QUIN & EP
            dSpine = 0.25 * (atp_pyr * (0.65 + 0.65 * sirt1_stac) + 1.15 * ketamine_burst) * (1.42 - S_spine) - 0.18 * (0.35 + 0.25 * QUIN + 0.012 * max(0.0, EP - 12.0)) * S_spine
            S_spine = max(0.45, min(1.38, S_spine + dt * dSpine))
            if abs(t - 7.0) < dt / 2:
                spine_day7 = round(S_spine, 2)

        return {
            "condition": condition,
            "arm": arm,
            "h_day90": round(D / (W + D) * 100, 1),
            "QUIN_day90": round(QUIN, 2),
            "GDF15_day90": round(GDF15, 0),
            "Hedonic_Drive_day90": hedonic_drive,
            "Spine_Density_day7": spine_day7,
            "Spine_Density_day90": round(S_spine, 2),
            "ATP_day90": round(atp_pyr, 2),
            "switch_day": switch_day
        }

    elif condition == "Parkinsons_SNc":
        W, D = 0.42, 0.58  # 58% baseline SNc mtDNA deletion / Complex I damaged fraction
        H2S = 0.58         # Desulfovibrio H2S + Enterobacteriaceae Curli vagal input
        DA_Q = 1.45        # Cytosolic Dopamine-o-Quinone (VMAT2 ATP deficit + Fe2+/O2 oxidation)
        Auto_Jam = 0.55    # GBA1 Asn370Ser + LRRK2 Gly2019Ser autophagosome-lysosome traffic jam
        SNCA_agg = 1.38    # Alpha-synuclein protofibril load on TOMM20/TOMM40
        pO2_snc = 1.72     # Mootha/Ichinose Brain Tissue Hyperoxia (1.00 = normoxic; >1.40 = O2 mismatch)
        GDF15 = 1950.0
        f_snc = 1.85       # Autonomous CaV1.3 pacemaker flux (f^2 = 3.42x)
        tissue_scale = 0.62
        stage = 1
        switch_day = None
        jam_day25 = None

        for step in range(steps):
            t = step * dt
            h_frac = D / max(W + D, 1e-6)
            if arm in ("PD_BCR_Nigrostriatal_Protocol", "PD_BCR_Plus_GBT601_PT2399_Shield"):
                gate_max_day = 26.5 if arm == "PD_BCR_Plus_GBT601_PT2399_Shield" else 34.0
                if stage == 1 and t > 10.0 and (h_frac <= 0.085 or GDF15 < 460.0 or t >= gate_max_day):
                    stage = 2
                    switch_day = round(t, 1)

            gbt601_pt2399_clamp = 0.0
            if arm == "Standard_Levodopa_Care_Alone":
                u_uro, u_nad, gba1_ambroxol, dbhb, tomm40_open, sirt1_stac = 0.08, 0.12, 0.10, 0.0, 0.22, 0.10
                ldopa_push = 0.42  # Unbuffered cytosolic dopamine load when VMAT2 lacks ATP
                dH2S = 0.16 - 0.28 * H2S
                rep_adv = 1.32
            elif arm == "Mitophagy_Alone_Without_GBA1_Ambroxol":
                # Discovery 6 Control: Stimulating PINK1/Parkin mitophagy without Ambroxol/ATP13A2 lysosomal rescue
                u_uro = 0.85
                u_nad, gba1_ambroxol, dbhb, tomm40_open, sirt1_stac = 0.22, 0.12, 0.15, 0.32, 0.20
                ldopa_push = 0.25
                dH2S = 0.14 - 0.30 * H2S
                rep_adv = 1.26
            elif arm == "Hypoxia_Pill_GBT601_PT2399_Alone":
                # Discovery 9 (Mootha & Ichinose Nature Neurosci / PNAS): GBT601 (Hb-O2 clamp, PDB: 5E83)
                # + PT2399 (HIF-2a/EPAS1 PAS-B inhibitor, PDB: 5UFP) monotherapy without GBA1/mitophagy clearance.
                u_uro, u_nad, gba1_ambroxol, dbhb, tomm40_open, sirt1_stac = 0.10, 0.14, 0.10, 0.15, 0.30, 0.14
                gbt601_pt2399_clamp = 0.92
                ldopa_push = 0.22
                dH2S = 0.15 - 0.28 * H2S
                rep_adv = 1.26
            elif arm in ("PD_BCR_Nigrostriatal_Protocol", "PD_BCR_Plus_GBT601_PT2399_Shield"):
                if stage == 1:
                    in_fast = ((t * 24.0) % 24.0) < 18.0
                    eff_p = max(0.0, ((t * 24.0) % 24.0) - 4.0) / 14.0 if in_fast else 0.0
                    u_uro = 0.96 if eff_p > 0 else 0.32
                    u_nad = 0.12
                    sirt1_stac = 0.25
                    rep_adv = 1.18
                    gbt601_pt2399_clamp = 0.94 if arm == "PD_BCR_Plus_GBT601_PT2399_Shield" else 0.0
                else:
                    in_p = (t % 7.0) < 1.75
                    u_uro = 0.90 if in_p else 0.22
                    u_nad = 0.15 if in_p else 0.98  # NADPARK 1000mg NR + PncA Synbiotic
                    sirt1_stac = 0.95               # 4ZZJ SIRT1 STAC + CoQ10 + Creatine
                    rep_adv = 1.04
                    gbt601_pt2399_clamp = 0.45 if arm == "PD_BCR_Plus_GBT601_PT2399_Shield" else 0.0
                gba1_ambroxol = 0.94  # 2NSX Ambroxol chaperone + Spermidine (ATP13A2) + Trehalose + TUDCA
                dbhb = 0.88           # C8-MCT + D-BHB Complex II bypass & VMAT2 ATP rescue
                tomm40_open = 0.92    # 7CK6 Citicoline + LPC-DHA + Sulforaphane DJ-1 rescue
                ldopa_push = 0.08     # Protected by Enterococcus TyrDC inhibition + VMAT2 sequestration
                dH2S = 0.035 - 0.38 * H2S

            H2S = max(0.05, H2S + dt * dH2S)

            # Autophagosome-lysosome traffic jam (Auto_Jam):
            # Increases when mitophagy (u_uro) pushes cargo into GBA1/LRRK2-blocked lysosomes; drained by Ambroxol + Spermidine (gba1_ambroxol)
            dJam = (0.34 * u_uro * (1.0 - 0.78 * gba1_ambroxol) + 0.14 * SNCA_agg * (1.0 - 0.65 * gba1_ambroxol)) - (0.28 + 0.72 * gba1_ambroxol) * Auto_Jam
            Auto_Jam = max(0.08, min(1.85, Auto_Jam + dt * dJam))
            if abs(t - 25.0) < dt / 2:
                jam_day25 = round(Auto_Jam, 2)

            # Effective lysosomal flux: choked if Auto_Jam is high
            lyso_flux = max(0.08, gba1_ambroxol / (1.0 + 1.35 * max(0.0, Auto_Jam - 0.30)))

            # OXPHOS & VMAT2 vesicular ATP
            eff_oxphos = max(0.16, (W + 0.12 * D) * (0.58 + 0.42 * tomm40_open) * max(0.35, 1.0 - 0.22 * DA_Q))
            # Tissue O2 Mismatch / Brain Hyperoxia (Mootha & Ichinose): rises when eff_oxphos < 0.85, clamped by GBT601+PT2399
            target_po2 = max(0.92, (1.0 + 0.95 * max(0.0, 0.88 - eff_oxphos)) * (1.0 - 0.44 * gbt601_pt2399_clamp))
            pO2_snc = max(0.90, min(2.10, pO2_snc + 0.35 * dt * (target_po2 - pO2_snc)))

            atp_snc = eff_oxphos * (1.0 + 0.35 * dbhb + 0.16 * sirt1_stac + 0.22 * gbt601_pt2399_clamp) / (1.0 + 0.55 * (H2S ** 1.4))

            # Cytosolic Dopamine-o-Quinone (DA_Q): generated when VMAT2 lacks ATP or excess tissue O2 auto-oxidizes dopamine;
            # quenched when GBT601+PT2399 clamps hyperoxia or VMAT2 + DJ-1/Sulforaphane sequesters DA
            vmat2_sequestration = min(1.0, 0.30 + 0.60 * atp_snc + 0.22 * dbhb + 0.28 * gbt601_pt2399_clamp)
            hyperoxia_o2_accel = 1.0 - 0.68 * gbt601_pt2399_clamp
            dDAQ = (
                (0.45 * (1.18 - vmat2_sequestration) + 0.32 * ldopa_push + 0.24 * max(0.0, Auto_Jam - 0.35)) * hyperoxia_o2_accel
                - (0.35 + 0.42 * gba1_ambroxol + 0.28 * dbhb + 0.65 * gbt601_pt2399_clamp) * DA_Q
            )
            DA_Q = max(0.08, min(2.40, DA_Q + dt * dDAQ))

            # Alpha-synuclein protofibril load (SNCA_agg): seeded by vagal H2S/Curli + DA_Q; cleared by unblocked GBA1/ATP13A2 lysosomes
            hypoxia_only_snca_persist = 0.24 if arm == "Hypoxia_Pill_GBT601_PT2399_Alone" else 0.0
            dSNCA = (0.24 * H2S + 0.22 * DA_Q + 0.20 * max(0.0, Auto_Jam - 0.40) + hypoxia_only_snca_persist) - (0.32 + 0.58 * lyso_flux) * SNCA_agg
            SNCA_agg = max(0.12 if gbt601_pt2399_clamp > 0 else 0.15, min(2.40, SNCA_agg + dt * dSNCA))

            eR_gut = 1.0 + 0.85 * (H2S ** 1.45) + 0.40 * max(0.0, DA_Q - 0.35) + 0.30 * max(0.0, SNCA_agg - 0.45) * (1.0 - 0.55 * gbt601_pt2399_clamp)
            eR = 10.0 * ((1.0 / max(0.18, eff_oxphos)) ** 1.22) * eR_gut
            f_eff = f_snc * (1.0 - 0.20 * dbhb - 0.12 * gbt601_pt2399_clamp)  # D-BHB & normoxic pO2 calm CaV1.3 spiking
            EP = eR * (f_eff ** 2)
            GDF15 = max(180.0, GDF15 + dt * (230.0 + 26.0 * max(0.0, EP - 14.0) - 0.85 * GDF15))

            # Net mitophagy clearance requires BOTH PINK1/Parkin initiation (u_uro) AND unblocked GBA1/ATP13A2 lysosomal fusion (lyso_flux)
            clearance = (0.072 * u_uro * (0.12 + 0.95 * lyso_flux) * (1.0 + 0.18 * gbt601_pt2399_clamp)) * tissue_scale
            biogenesis = (0.024 + 0.165 * u_nad * (0.85 + 0.45 * sirt1_stac) * tomm40_open * max(0.35, 1.0 - 0.28 * SNCA_agg)) * tissue_scale
            crowding = max(0.0, 1.0 - (W + D) / 1.40)
            W = max(0.05, W + dt * (biogenesis * W * crowding - 0.008 * tissue_scale * W))
            D = max(0.001, D + dt * (biogenesis * rep_adv * D * crowding - (0.008 * tissue_scale + clearance * (1.55 + 1.15 * (D / (W + D)))) * D))

        # Functional Nigrostriatal Vesicular Dopamine Tone (1.00 = healthy SNc striatal release)
        # Note (Mootha & Ichinose): When GBT601+PT2399 clamps hyperoxia & DA_Q, surviving SNc neurons restore dopamine tone even if SNCA_agg remains!
        snca_penalty = 0.35 * max(0.0, SNCA_agg - 0.35) * (1.0 - 0.80 * gbt601_pt2399_clamp)
        max_th = 1.28 if gbt601_pt2399_clamp > 0 else 1.25
        max_res = 98.2 if gbt601_pt2399_clamp > 0 else 95.0
        th_tone = round(min(max_th, (0.28 + 0.65 * atp_snc + 0.36 * gbt601_pt2399_clamp) * vmat2_sequestration / (1.0 + snca_penalty)), 2)
        updrs_improvement = round(max(0.0, min(max_res, (th_tone - 0.24) / (1.05 - 0.24) * 86.0)), 1)
        return {
            "condition": condition,
            "arm": arm,
            "h_day90": round(D / (W + D) * 100, 1),
            "H2S_day90": round(H2S, 2),
            "Auto_Jam_day25": jam_day25,
            "Auto_Jam_day90": round(Auto_Jam, 2),
            "pO2_SNc_day90": round(pO2_snc, 2),
            "DA_Q_day90": round(DA_Q, 2),
            "SNCA_agg_day90": round(SNCA_agg, 2),
            "GDF15_day90": round(GDF15, 0),
            "ATP_SNc_day90": round(atp_snc, 2),
            "TH_Dopamine_Tone_day90": th_tone,
            "Motor_Score_Rescue_Pct": updrs_improvement,
            "switch_day": switch_day
        }

print("\n=== DISEASE-SPECIFIC SUB-PROTOCOL DEEP DIVE: EPILEPSY, TR-MDD & PARKINSON'S DISEASE (SNc) ===")
neuro_results = {"epilepsy": [], "mdd": [], "parkinsons": []}
for arm in ["Standard_Anticonvulsant_Care", "Chronic_Classic_Ketogenic_Diet", "Epi_BCR_Protocol"]:
    r = simulate_neuropsychiatric_subprotocols("Refractory_Epilepsy", arm)
    neuro_results["epilepsy"].append(r)
    print(f"  [Epilepsy]   {arm:38s} | h={r['h_day90']:4.1f}% | H2S={r['H2S_day90']:.2f} | PV+ ATP={r['ATP_PV_day90']:.2f}x | E/I Ratio d25={r['EI_ratio_day25']:.2f} -> d90={r['EI_ratio_day90']:.2f} | Seizure Red={r['Seizure_Reduction_Pct']:4.1f}% (Sw d{r['switch_day']})")

for arm in ["Standard_SSRI_Monoamine_Push", "Acute_Ketamine_Burst_Alone", "Affective_MDD_BCR_Protocol"]:
    r = simulate_neuropsychiatric_subprotocols("Treatment_Resistant_MDD", arm)
    neuro_results["mdd"].append(r)
    print(f"  [TR-MDD]     {arm:38s} | h={r['h_day90']:4.1f}% | QUIN={r['QUIN_day90']:.2f} | GDF15={r['GDF15_day90']:4.0f} | Drive={r['Hedonic_Drive_day90']:.2f}x | Spines d7={r['Spine_Density_day7']:.2f}x -> d90={r['Spine_Density_day90']:.2f}x (Sw d{r['switch_day']})")

for arm in [
    "Standard_Levodopa_Care_Alone",
    "Mitophagy_Alone_Without_GBA1_Ambroxol",
    "Hypoxia_Pill_GBT601_PT2399_Alone",
    "PD_BCR_Nigrostriatal_Protocol",
    "PD_BCR_Plus_GBT601_PT2399_Shield",
]:
    r = simulate_neuropsychiatric_subprotocols("Parkinsons_SNc", arm)
    neuro_results["parkinsons"].append(r)
    print(f"  [Parkinsons] {arm:38s} | h={r['h_day90']:4.1f}% | pO2={r['pO2_SNc_day90']:.2f}x | Jam d25={r['Auto_Jam_day25']:.2f}->d90={r['Auto_Jam_day90']:.2f} | DA-Q={r['DA_Q_day90']:.2f} | SNCA={r['SNCA_agg_day90']:.2f} | TH Tone={r['TH_Dopamine_Tone_day90']:.2f}x | Rescue={r['Motor_Score_Rescue_Pct']:4.1f}% (Sw d{r['switch_day']})")

with open("/Users/o2satz/.gemini/antigravity/brain/e55fbd0d-7b92-471e-b9ba-63a180b98a1f/scratch/data/bcr_neuropsychiatric_results.json", "w") as f:
    json.dump(neuro_results, f, indent=2)


def simulate_cart_bioreactor(arm="3_Bioreactor_BCR_Biphasic_Gated_Controller"):
    """
    Discovery 7 (Bioreactor-BCR): 10-Day (240-Hour) Ex Vivo CAR-T Cell Therapy Manufacturing Controller.
    Models:
      - Wild-type (W) vs. deletion/depolarized (D) mitochondrial mass per T-cell during clonal expansion
      - Pyruvate Dehydrogenase Kinase 1 (PDK1, PDB: 2Q8G) Warburg lock vs. AZD7545/DCA L2-pocket inhibition (Ser75 2.50 A)
      - Automated Hour-60 Supernatant Sensor Gate (Lactate:Pyruvate < 12:1 & GDF15 < 380 pg/mL)
      - Stem-like Memory T-cell fraction (T_SCM, TCF7+/LEF1+/CCR7+/CD62L+), Spare Respiratory Capacity (SRC), and Fold Expansion
    """
    dt_h = 0.1  # 6-minute resolution across 240 hours (10 days)
    steps = int(240.0 / dt_h)
    W, D = 0.48, 0.42          # Patient apheresis T-cells starting at 46.7% damaged mitochondria (post-chemo oncology patient)
    cell_fold = 1.0            # Normalized starting CAR-T cell count (1.0x)
    t_scm = 0.42               # Starting naive/stem-like TCF7+ fraction
    exhaustion = 0.28          # Starting TIM-3+/PD-1+/TOX+ exhaustion marker fraction
    lac_pyr = 22.0             # Supernatant Lactate:Pyruvate ratio
    gdf15 = 820.0              # Supernatant GDF15 (pg/mL)
    stage = 1
    switch_hour = None

    for step in range(steps):
        hr = step * dt_h
        dt_d = dt_h / 24.0
        h_frac = D / max(W + D, 1e-6)

        if arm == "1_Standard_CD3_CD28_IL2_HighGlucose":
            # Continuous 25 mM glucose + IL-2 + strong CD3/CD28: locks mTORC1 ON and PDK1 (2Q8G) active
            pdk1_inhib = 0.05
            mtorc1 = 0.94
            u_uro_ktp = 0.08
            u_nad_stac = 0.15
            il7_il15 = 0.10
            rep_adv = 1.38  # High Warburg reverse-ATPase masks D while POLG rapidly clones D
        elif arm == "2_Static_Continuous_UrolithinA_Rapamycin":
            # Continuous 10-day mTOR inhibition + Urolithin A without Phase 2 anabolic release
            pdk1_inhib = 0.80
            mtorc1 = 0.22
            u_uro_ktp = 0.90
            u_nad_stac = 0.30
            il7_il15 = 0.75
            rep_adv = 1.05
        elif arm == "3_Bioreactor_BCR_Biphasic_Gated_Controller":
            # Phase 1 (Hours 0-60): Low-glucose/galactose + AZD7545/DCA PDK1 block (2Q8G) + Urolithin A + Kinetin (6EQI)
            # Transitions to Phase 2 when h_frac <= 0.065 or (hr >= 54 and lac_pyr < 12.0 and gdf15 < 380.0)
            if stage == 1 and hr >= 48.0 and (h_frac <= 0.065 or (lac_pyr < 12.0 and gdf15 < 380.0) or hr >= 64.0):
                stage = 2
                switch_hour = round(hr, 1)
            if stage == 1:
                pdk1_inhib = 0.94  # 2Q8G Ser75/Phe62 lipoyl-pocket block unblocks PDH flux
                mtorc1 = 0.24
                u_uro_ktp = 0.96   # Urolithin A + Kinetin (6EQI) mitophagy pulse
                u_nad_stac = 0.20
                il7_il15 = 0.85
                rep_adv = 1.06
            else:
                pdk1_inhib = 0.82
                mtorc1 = 0.86      # Controlled anabolic perfusion release for clonal expansion
                u_uro_ktp = 0.22
                u_nad_stac = 0.96  # Nicotinamide Riboside + SIRT1 STAC (4ZZJ) + Citicoline (7CK6) + Creatine/Arginine
                il7_il15 = 0.95
                rep_adv = 1.02

        # Supernatant Lactate:Pyruvate ratio driven by Warburg glycolysis vs PDK1 inhibition
        target_lp = 28.5 * (1.0 - 0.68 * pdk1_inhib) + 2.0 * (1.0 - pdk1_inhib)
        lac_pyr = lac_pyr + (target_lp - lac_pyr) * 0.85 * dt_d

        # Voltage unmasking of damaged mitochondria requires low Lactate:Pyruvate (PDH active, reverse-ATPase off)
        voltage_unmask = 1.0 / (1.0 + math.exp((lac_pyr - 14.5) / 2.5))
        ulk1_open = 1.0 - 0.70 * mtorc1

        # Rapid T-cell mitochondrial turnover (tissue_scale = 2.6x in activated lymphocytes)
        t_scale = 2.6
        clearance_D = t_scale * (0.06 + 0.55 * u_uro_ktp * voltage_unmask * ulk1_open)
        bio_W = t_scale * (0.05 + 0.30 * (0.30 + 0.70 * mtorc1) * (0.45 + 0.80 * u_nad_stac))
        bio_D = t_scale * (0.05 + 0.30 * (0.30 + 0.70 * mtorc1) * (0.45 + 0.55 * u_nad_stac)) * rep_adv * (1.0 - 0.48 * voltage_unmask)

        # Reductive stress EP = eR * f^2 in bioreactor
        eff_oxphos = max(0.22, W + 0.18 * D)
        eR = 7.5 * ((1.0 / eff_oxphos) ** 1.15)
        f_push = 1.28 if arm == "1_Standard_CD3_CD28_IL2_HighGlucose" else (0.82 if stage == 1 else 1.08)
        EP = eR * (f_push ** 2)
        dmg = 0.008 * max(0.0, EP - 9.5) * W

        crowding = max(0.0, 1.0 - (W + D) / 1.48)
        W = max(0.12, W + dt_d * (bio_W * W * crowding - 0.02 * W - dmg))
        D = max(0.002, D + dt_d * (bio_D * D * crowding - (0.02 + clearance_D) * D + dmg))

        target_gdf15 = 185.0 + 22.0 * EP
        gdf15 = gdf15 + (target_gdf15 - gdf15) * 0.90 * dt_d

        # Clonal division rate: Warburg glycolysis (mtorc1) drives rapid division initially, while OXPHOS + mtorc1 in Phase 2 drives superior sustained expansion
        src_curr = (W + 0.15 * D) * (1.0 + 0.55 * pdk1_inhib + 0.65 * u_nad_stac) * 1.28
        if arm == "1_Standard_CD3_CD28_IL2_HighGlucose":
            div_rate = 0.395  # Standard glycolytic expansion (~52x over 10 days), but clones D!
        elif arm == "2_Static_Continuous_UrolithinA_Rapamycin":
            div_rate = 0.238  # Chronic mTORC1 blockade stunts total division (~10.8x over 10 days)
        else:
            div_rate = 0.255 if stage == 1 else 0.472  # Phase 1 purge followed by rapid Phase 2 OXPHOS-fueled division (~64x)
        cell_fold = cell_fold * math.exp(div_rate * dt_d)

        # Stem-like memory T_SCM (TCF7+/LEF1+) promoted by alpha-KG/OXPHOS + u_uro_ktp + IL-7/15; eroded by high Warburg EP
        target_tscm = min(0.86, max(0.10, 0.24 + 0.38 * pdk1_inhib + 0.26 * il7_il15 - 0.45 * h_frac))
        t_scm = t_scm + (target_tscm - t_scm) * 0.45 * dt_d

        # Exhaustion (TIM-3+ / PD-1+ / TOX+) driven by EP and heteroplasmy
        target_exh = min(0.85, max(0.04, 0.06 + 0.82 * h_frac + 0.012 * max(0.0, EP - 9.0) - 0.22 * u_nad_stac))
        exhaustion = exhaustion + (target_exh - exhaustion) * 0.45 * dt_d

    src_final = (W + 0.15 * D) * (1.0 + 0.55 * pdk1_inhib + 0.65 * u_nad_stac) * 1.28
    return {
        "arm": arm,
        "switch_hour": switch_hour,
        "h_hour240_pct": round(D / (W + D) * 100.0, 1),
        "W_hour240": round(W, 2),
        "Lactate_Pyruvate_hour240": round(lac_pyr, 1),
        "GDF15_hour240_pg_ml": round(gdf15, 0),
        "Fold_Expansion_day10_x": round(cell_fold, 1),
        "T_SCM_Stemness_pct": round(t_scm * 100.0, 1),
        "Exhaustion_TIM3_PD1_pct": round(exhaustion * 100.0, 1),
        "Spare_Respiratory_Capacity_x": round(src_final, 2)
    }


def simulate_apis_honeybee(arm="3_Apis_BCR_Two_Stage_Holobiont_Patty"):
    """
    Discovery 8 (Apis-BCR): 14-Day Honeybee (Apis mellifera) Flight-Muscle & Snodgrassella alvi Holobiont Rescue.
    Models:
      - Strobilurin Fungicide (Azoxystrobin / Pyraclostrobin, PDB: 1SQB) Qo-site blockade (Glu271 3.15 A, Phe274 3.33 A)
      - Neonicotinoid (Imidacloprid) Complex IV (cox5a/cox6c) suppression & Glyphosate Snodgrassella alvi dysbiosis
      - Wingbeat Transition at Day 6: In-Hive Nurse/Resting flux (f_rest = 1.0) -> Foraging Flight (100-250 Hz, f_flight = 2.40 => f^2 = 5.76x)
      - Stage 1 (Days 1-5): p-Coumaric Acid (CYP9Q1/2/3 detox) + Coumaroyl-Spermidine/Urolithin A mitophagy + Live S. alvi synbiotic
      - Stage 2 (Days 6-14): Ubiquinol/CoQ10 (competitive Qo-site rescue at Ala143/Ile146 in 1SQB) + NR + Riboflavin (B2)
    """
    dt = 0.02  # 28.8-minute resolution across 14 days
    steps = int(14.0 / dt)
    W, D = 0.48, 0.52       # 52% damaged thoracic flight-muscle mitochondria after agricultural fungicide + neonicotinoid exposure
    strob_bound = 0.78      # 78% of Complex III Qo sites (PDB: 1SQB) initially blocked by Azoxystrobin/Pyraclostrobin
    s_alvi = 0.28           # Depleted Snodgrassella alvi ileum symbiont density (normal = 1.0)
    stage = 1
    switch_day = None
    ep_flight_peak = 0.0

    for step in range(steps):
        t = step * dt
        h_frac = D / max(W + D, 1e-6)

        if arm == "1_Exposed_Plain_HFCS_Syrup":
            p_coumaric_cyp9q = 0.08
            spermidine_uroA = 0.06
            salvi_synbiotic = 0.05
            coq10_nr_b2 = 0.08
            hfcs_glycolytic_push = 1.15  # High-fructose corn syrup forces unbuffered glycolysis into blocked Complex III
            rep_adv = 1.32
        elif arm == "2_Unphased_Antioxidant_Pollen_Sub_Without_Detox_Purge":
            # CoQ10/vitamins fed continuously without p-Coumaric CYP9Q detox or Stage 1 Coumaroyl-Spermidine mitophagy purge
            p_coumaric_cyp9q = 0.14
            spermidine_uroA = 0.08
            salvi_synbiotic = 0.42
            coq10_nr_b2 = 0.55
            hfcs_glycolytic_push = 1.05
            rep_adv = 1.28
        elif arm == "3_Apis_BCR_Two_Stage_Holobiont_Patty":
            if stage == 1 and t >= 4.5 and (h_frac <= 0.075 or strob_bound <= 0.18 or t >= 5.2):
                stage = 2
                switch_day = round(t, 1)
            if stage == 1:
                # Stage 1 (Days 1-5 In-Hive Nurse/Resting Detox & Mitophagy Patty)
                p_coumaric_cyp9q = 0.96   # Floral p-coumaric acid upregulates CYP9Q1/2/3 to clear strobilurin/neonicotinoid
                spermidine_uroA = 0.96    # Coumaroyl-spermidine + Urolithin A purges damaged flight-muscle mitochondria
                salvi_synbiotic = 0.95    # Live Snodgrassella alvi + Gilliamella apicola + organic acid prebiotic
                coq10_nr_b2 = 0.25
                hfcs_glycolytic_push = 0.72
                rep_adv = 1.06
            else:
                # Stage 2 (Days 6-14 Forager Flight OXPHOS Build Patty)
                p_coumaric_cyp9q = 0.88
                spermidine_uroA = 0.28
                salvi_synbiotic = 0.95
                coq10_nr_b2 = 0.98        # Ubiquinol (CoQ10) + NR + Riboflavin B2 + 24-Methylenecholesterol
                hfcs_glycolytic_push = 0.92
                rep_adv = 1.02

        # 1. CYP9Q1/2/3 P450 clearance of Strobilurin from Complex III Qo site (PDB: 1SQB) + CoQ10 competitive displacement
        dStrob = 0.03 - (0.38 * p_coumaric_cyp9q + 0.18 * coq10_nr_b2) * strob_bound
        strob_bound = max(0.04, min(0.95, strob_bound + dt * dStrob))

        # 2. Snodgrassella alvi ileum colonization recovery (restores gut O2 barrier & citrate/malate pool)
        dSalvi = (0.35 * salvi_synbiotic * (1.02 - s_alvi)) - (0.12 * (1.0 - salvi_synbiotic) * s_alvi)
        s_alvi = max(0.15, min(1.00, s_alvi + dt * dSalvi))

        # 3. Thoracic Flight-Muscle Energy Resistance (eR) & Foraging Wingbeat Flux (f_flight = 2.40 => f^2 = 5.76x after Day 5.5)
        in_foraging_flight = (t >= 5.5)
        f_wing = (2.40 if in_foraging_flight else 1.00) * hfcs_glycolytic_push
        q_site_conductance = max(0.28, (1.0 - 0.68 * strob_bound) * (0.65 + 0.45 * coq10_nr_b2))
        gut_factor = 1.0 + 0.45 * (1.0 - s_alvi)
        eff_oxphos = max(0.22, (W + 0.15 * D) * q_site_conductance)
        eR = 2.8 * ((1.0 / eff_oxphos) ** 1.05) * gut_factor
        EP = eR * (f_wing ** 2)
        if in_foraging_flight and EP > ep_flight_peak:
            ep_flight_peak = round(EP, 1)

        # 4. Flight-muscle mitophagy & biogenesis (faster insect metabolic rate, scale = 2.2x)
        bee_scale = 2.2
        clearance_D = bee_scale * (0.04 + 0.38 * spermidine_uroA)
        bio_W = bee_scale * (0.035 + 0.22 * coq10_nr_b2 * s_alvi)
        bio_D = bee_scale * (0.035 + 0.16 * coq10_nr_b2) * rep_adv * (1.0 - 0.48 * spermidine_uroA)
        dmg = 0.0022 * max(0.0, EP - 16.0) * W

        crowding = max(0.0, 1.0 - (W + D) / 1.42)
        W = max(0.14, W + dt * (bio_W * W * crowding - 0.015 * W - dmg))
        D = max(0.002, D + dt * (bio_D * D * crowding - (0.015 + clearance_D) * D + dmg))

    flight_atp = round((W + 0.15 * D) * q_site_conductance * (1.0 + 0.22 * s_alvi) / gut_factor, 2)
    homing_survival = round(max(18.0, min(96.5, 24.0 + 58.0 * min(1.3, flight_atp) - 0.32 * max(0.0, ep_flight_peak - 18.0))), 1)
    return {
        "arm": arm,
        "switch_day": switch_day,
        "h_day14_pct": round(D / (W + D) * 100.0, 1),
        "Strobilurin_Qo_Block_day14_pct": round(strob_bound * 100.0, 1),
        "Snodgrassella_alvi_Colonization_x": round(s_alvi, 2),
        "Peak_Flight_EP_Stress": ep_flight_peak,
        "Thoracic_Flight_ATP_day14_x": flight_atp,
        "Forager_Homing_Survival_pct": homing_survival
    }


print("\n=== DISCOVERY 7: EX VIVO CAR-T MANUFACTURING CONTROLLER (Cohort F: Bioreactor-BCR | 240h / 10 Days | PDB: 2Q8G) ===")
cart_results = []
for arm in [
    "1_Standard_CD3_CD28_IL2_HighGlucose",
    "2_Static_Continuous_UrolithinA_Rapamycin",
    "3_Bioreactor_BCR_Biphasic_Gated_Controller"
]:
    r = simulate_cart_bioreactor(arm)
    cart_results.append(r)
    print(f"  [CAR-T]      {arm:44s} | h={r['h_hour240_pct']:4.1f}% | Lac:Pyr={r['Lactate_Pyruvate_hour240']:4.1f} | Fold={r['Fold_Expansion_day10_x']:4.1f}x | T_SCM={r['T_SCM_Stemness_pct']:4.1f}% | SRC={r['Spare_Respiratory_Capacity_x']:.2f}x | Exh={r['Exhaustion_TIM3_PD1_pct']:4.1f}% (Sw h{r['switch_hour']})")

print("\n=== DISCOVERY 8: HONEYBEE (APIS MELLIFERA) FLIGHT-MUSCLE & HOLOBIONT RESCUE (Cohort G: Apis-BCR | 14 Days | PDB: 1SQB) ===")
apis_results = []
for arm in [
    "1_Exposed_Plain_HFCS_Syrup",
    "2_Unphased_Antioxidant_Pollen_Sub_Without_Detox_Purge",
    "3_Apis_BCR_Two_Stage_Holobiont_Patty"
]:
    r = simulate_apis_honeybee(arm)
    apis_results.append(r)
    print(f"  [Apis-BCR]   {arm:52s} | h={r['h_day14_pct']:4.1f}% | Qo-Block={r['Strobilurin_Qo_Block_day14_pct']:4.1f}% | S.alvi={r['Snodgrassella_alvi_Colonization_x']:.2f}x | Flight EP={r['Peak_Flight_EP_Stress']:4.1f} | ATP={r['Thoracic_Flight_ATP_day14_x']:.2f}x | Homing={r['Forager_Homing_Survival_pct']:4.1f}% (Sw d{r['switch_day']})")

with open("/Users/o2satz/.gemini/antigravity/scratch/mitorestorer/data/bcr_cross_species_cart_apis.json", "w") as f:
    json.dump({"cart_bioreactor": cart_results, "apis_honeybee": apis_results}, f, indent=2)





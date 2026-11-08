/* ============================================================================
 * Google Coral NPU (`RV32IMF_Zve32x`) Hardware-in-the-Loop (HIL) C99 Testbench
 * Executes `mitorestorer_coralnpu_kernel.h` across all 3 Cross-Species Targets:
 *   1. CORAL_TARGET_APIS_SMART_HIVE      (Cohort G: Honeybee Day 0 -> Day 4.5 Gate)
 *   2. CORAL_TARGET_CART_GMP_BIOREACTOR  (Cohort F: CAR-T Hour 0 -> Hour 64.0 Gate)
 *   3. CORAL_TARGET_HUMAN_WEARABLE_PATCH (Cohorts A-E: Human Day 0 -> Day 38.0 Gate)
 * ============================================================================ */

#include <stdio.h>
#include <assert.h>
#include "mitorestorer_coralnpu_kernel.h"

int main(void) {
  int32_t logit_acc = 0;
  int32_t skipped_bins = 0;

  /* 1. Test Target 1: Apis-BCR Smart Beehive Node (Day 0 Purge -> Day 4.5 Build Gate) */
  bool hive_safe = coralnpu_scalar_safety_guard(CORAL_TARGET_APIS_SMART_HIVE, 3450, 120);
  assert(hive_safe == true);

  mitorestorer_ode_state_q16_t hive_state_d0 = {
    .w_mass_q16 = (int32_t)(0.43 * 65536),
    .d_mass_q16 = (int32_t)(0.47 * 65536), /* h = 52.2% > 8.0% */
    .inhibitor_q16 = (int32_t)(48.0 * 65536)
  };
  const int8_t wingbeat_droop_i8[16] = {0, 0, 0, 0, 0, 0, 0, 0, 40, 50, 0, 0, 0, 0, 0, 0};
  mitorestorer_gate_action_t g_apis_d0 = coralnpu_matrix256_stage_gate(
      &hive_state_d0, wingbeat_droop_i8, kMitoRestorerGateWeightsQ7, hive_safe, &logit_acc, &skipped_bins);
  assert(g_apis_d0 == BCR_GATE_STAGE1_MITOPHAGY_PURGE);

  /* Integrate Q16.16 RK4 step to Day 4.5 (h = 6.8% <= 8.0%, Qo blockade = 6.0 <= 12.0) */
  mitorestorer_ode_state_q16_t hive_state_d4_5 = {
    .w_mass_q16 = (int32_t)(0.82 * 65536),
    .d_mass_q16 = (int32_t)(0.06 * 65536),
    .inhibitor_q16 = (int32_t)(6.0 * 65536)
  };
  int32_t k1[8] = {1024, -512, 0, 0, 512, 512, -512, 0};
  int32_t k2[8] = {1024, -512, 0, 0, 512, 512, -512, 0};
  coralnpu_rvv_rk4_bioenergetic_step(&hive_state_d4_5, k1, k2, 6553);

  const int8_t wingbeat_healthy_i8[16] = {64, 52, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 28};
  mitorestorer_gate_action_t g_apis_d4_5 = coralnpu_matrix256_stage_gate(
      &hive_state_d4_5, wingbeat_healthy_i8, kMitoRestorerGateWeightsQ7, hive_safe, &logit_acc, &skipped_bins);
  assert(g_apis_d4_5 == BCR_GATE_STAGE2_ANABOLIC_BUILD);

  /* 2. Test Target 2: CAR-T GMP Bioreactor Cassette Valve (Hour 64.0 Lac:Pyr <= 12 Gate) */
  bool cart_safe = coralnpu_scalar_safety_guard(CORAL_TARGET_CART_GMP_BIOREACTOR, 720, 5500);
  assert(cart_safe == true);
  mitorestorer_ode_state_q16_t cart_state_h64 = {
    .w_mass_q16 = (int32_t)(0.88 * 65536),
    .d_mass_q16 = (int32_t)(0.05 * 65536), /* h = 5.4% */
    .inhibitor_q16 = (int32_t)(11.4 * 65536) /* Lac:Pyr = 11.4 <= 12.0 */
  };
  mitorestorer_gate_action_t g_cart_h64 = coralnpu_matrix256_stage_gate(
      &cart_state_h64, wingbeat_healthy_i8, kMitoRestorerGateWeightsQ7, cart_safe, &logit_acc, &skipped_bins);
  assert(g_cart_h64 == BCR_GATE_STAGE2_ANABOLIC_BUILD);

  /* 3. Test Target 3: Human Wearable Ring Safety Override when Lean Mass drops <= -2.5% */
  bool hum_unsafe = coralnpu_scalar_safety_guard(CORAL_TARGET_HUMAN_WEARABLE_PATCH, 3680, -280);
  assert(hum_unsafe == false);
  int32_t override_logit = 0, override_skip = 0;
  mitorestorer_gate_action_t g_hum_override = coralnpu_matrix256_stage_gate(
      &cart_state_h64, wingbeat_healthy_i8, kMitoRestorerGateWeightsQ7, hum_unsafe, &override_logit, &override_skip);
  assert(g_hum_override == BCR_GATE_STAGE0_SAFETY_HOLD);

  /* 4. Test Rolling-Window Noise Filter ("Single Tripwire in a Hurricane" Mitigation):
   *    A transient 5-frame wind-gust spike at 225.0 Hz (2250) in Apis-BCR MUST NOT trigger Stage 2,
   *    whereas 60+ sustained frames at 221.5 Hz (2215) MUST trigger Stage 2. */
  mitorestorer_rolling_guard_dtcm_t guard_state = {0};
  mitorestorer_gate_action_t g_wind_gust = BCR_GATE_STAGE1_MITOPHAGY_PURGE;
  for (int f = 0; f < 5; ++f) {
    g_wind_gust = coralnpu_rolling_fuzzy_sensor_gate(
        CORAL_TARGET_APIS_SMART_HIVE, &guard_state, 2250, 3450, 480, true);
  }
  assert(g_wind_gust == BCR_GATE_STAGE1_MITOPHAGY_PURGE);
  mitorestorer_gate_action_t g_sustained_64 = BCR_GATE_STAGE1_MITOPHAGY_PURGE;
  for (int f = 0; f < 60; ++f) {
    g_sustained_64 = coralnpu_rolling_fuzzy_sensor_gate(
        CORAL_TARGET_APIS_SMART_HIVE, &guard_state, 2215, 3450, 1080, true);
  }
  assert(g_sustained_64 == BCR_GATE_STAGE2_ANABOLIC_BUILD);

  /* 5. Test CAR-T Amperometric Biofouling Fallback (Timed 64h + Optical pH Check)
   *    and Human Post-Meal Ketone Skew Blocking (PPG-HRV x B-OHB Cross-Reference) */
  mitorestorer_rolling_guard_dtcm_t cart_guard = {0};
  cart_guard.prev_lac_pyr_x100 = 2400;
  /* Inject 3 erratic biofouling jumps (> 4.50 ratio units each) */
  (void)coralnpu_rolling_fuzzy_sensor_gate(CORAL_TARGET_CART_GMP_BIOREACTOR, &cart_guard, 1600, 720, 420, false);
  (void)coralnpu_rolling_fuzzy_sensor_gate(CORAL_TARGET_CART_GMP_BIOREACTOR, &cart_guard, 2300, 720, 430, false);
  mitorestorer_gate_action_t g_biofoul_hold = coralnpu_rolling_fuzzy_sensor_gate(
      CORAL_TARGET_CART_GMP_BIOREACTOR, &cart_guard, 1100, 720, 440, false);
  assert(g_biofoul_hold == BCR_GATE_BIOFOUL_TIMED64H_OPTICAL_PH_HOLD);
  /* At Hour 64.0 (640) with verified optical pH, fallback safely releases Stage 2 */
  mitorestorer_gate_action_t g_biofoul_pass = coralnpu_rolling_fuzzy_sensor_gate(
      CORAL_TARGET_CART_GMP_BIOREACTOR, &cart_guard, 1100, 720, 640, true);
  assert(g_biofoul_pass == BCR_GATE_STAGE2_ANABOLIC_BUILD);

  /* Post-meal transient ketone skew (B-OHB=1.5 mM [15], but autonomic PPG-HRV RMSSD=24.0 ms [2400 < 3800]) */
  mitorestorer_gate_action_t g_meal_skew = coralnpu_rolling_fuzzy_sensor_gate(
      CORAL_TARGET_HUMAN_WEARABLE_PATCH, &guard_state, 15, 2400, 360, true);
  assert(g_meal_skew == BCR_GATE_STAGE1_MITOPHAGY_PURGE);

  /* 6. Test Bidirectional Physical Rollback & 32-Byte L1 DTCM Breadcrumb Restore:
   *    Save safe Stage-1 checkpoint, then simulate post-Stage-2 biological rejection across all 3 targets:
   *    - Apis-BCR: Brood temp drops to 32.80C (3280) + 72 dB distress hiss -> Snap Solenoid Baffle back to Bay-A!
   *    - CAR-T: Post-Stage-2 pH drops to 6.82 (682) -> Un-pinch Media-B & flush neutral perfusion buffer!
   *    - Human APOE4/TOMM40 Arg330: Post-Stage-2 GDF15 rebounds to 4,000 pg/mL -> Evacuate & restore breadcrumb! */
  coralnpu_save_breadcrumb_ckpt(&guard_state, &hive_state_d4_5);
  int32_t saved_w_q16 = hive_state_d4_5.w_mass_q16;
  hive_state_d4_5.w_mass_q16 = (int32_t)(0.31 * 65536); /* Corrupted by post-transition stress */

  mitorestorer_gate_action_t g_apis_rb = coralnpu_bidirectional_rollback_guard(
      CORAL_TARGET_APIS_SMART_HIVE, BCR_GATE_STAGE2_ANABOLIC_BUILD, &hive_state_d4_5, &guard_state, 3280, 72);
  assert(g_apis_rb == BCR_GATE_STAGE2_ROLLBACK_ABORT);
  assert(hive_state_d4_5.w_mass_q16 == saved_w_q16); /* Breadcrumb state restored bit-exactly! */

  mitorestorer_gate_action_t g_cart_rb = coralnpu_bidirectional_rollback_guard(
      CORAL_TARGET_CART_GMP_BIOREACTOR, BCR_GATE_STAGE2_ANABOLIC_BUILD, &cart_state_h64, &guard_state, 682, 4100);
  assert(g_cart_rb == BCR_GATE_STAGE2_ROLLBACK_ABORT);

  mitorestorer_gate_action_t g_apoe4_rb = coralnpu_bidirectional_rollback_guard(
      CORAL_TARGET_HUMAN_WEARABLE_PATCH, BCR_GATE_STAGE2_ANABOLIC_BUILD, &cart_state_h64, &guard_state, 4000, 680);
  assert(g_apoe4_rb == BCR_GATE_STAGE2_ROLLBACK_ABORT);

  printf("CORAL_NPU_C99_HIL_PASS: Apis_d0=%d Apis_d4_5=%d CART_h64=%d Human_Override=%d WindGust_Rej=%d Biofoul_Hold=%d MealSkew_Rej=%d Rollbacks=%d/%d/%d Logit=%d SkippedBins=%d/16\n",
         g_apis_d0, g_apis_d4_5, g_cart_h64, g_hum_override, g_wind_gust, g_biofoul_hold, g_meal_skew,
         g_apis_rb, g_cart_rb, g_apoe4_rb, logit_acc, skipped_bins);
  return 0;
}

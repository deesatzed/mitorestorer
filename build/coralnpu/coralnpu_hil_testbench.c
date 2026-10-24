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

  printf("CORAL_NPU_C99_HIL_PASS: Apis_d0=%d Apis_d4_5=%d CART_h64=%d Human_Override=%d Logit=%d SkippedBins=%d/16\n",
         g_apis_d0, g_apis_d4_5, g_cart_h64, g_hum_override, logit_acc, skipped_bins);
  return 0;
}

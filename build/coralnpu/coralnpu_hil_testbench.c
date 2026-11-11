/* ============================================================================
 * Google Coral NPU (`RV32IMF_Zve32x`) Host-Emulated C99 Verification Testbench
 * File: build/coralnpu/coralnpu_hil_testbench.c
 * Executes:
 *   1. Cross-target decision & safety guardrail assertions (Targets 1-3)
 *   2. 64-frame persistence + secondary-modality enforcement & biofouling hold
 *   3. Bidirectional actuator rollback & 32-byte reference anchor verification
 *   4. Full 1,000-step 8-state 4-stage RK4 Q16.16 vs. IEEE-754 FP64 parity audit
 * ============================================================================ */

#include <stdio.h>
#include <math.h>
#include <assert.h>
#include "mitorestorer_coralnpu_kernel.h"

/* IEEE-754 FP64 Reference 8-State ODE RHS Evaluator (Exact Match to Eqs. 4-11) */
static void eval_core8_rhs_fp64(
    const double x[8],
    double u0,
    double u1,
    double u2,
    double dxdt[8]) {
  const double K_cap = 1.250;
  const double k_rep_adv = 1.320;
  const double k_bio0 = 0.075;
  const double k_biomax = 0.140;
  const double k_mito_w = 0.015;
  const double k_mito_d_max = 0.115;
  const double k_dmg = 0.040;
  const double k_acc = 0.080;
  const double k_clr = 0.180;
  const double K_B = 0.250;
  const double k_E = 0.250;
  const double k_g_prod = 0.080;
  const double k_g_clr = 0.160;
  const double k_n_syn0 = 0.060;
  const double k_n_syn2 = 0.080;
  const double k_parp = 0.060;
  const double k_n_use = 0.095;
  const double k_oct0 = 0.040;
  const double k_oct2 = 0.100;
  const double k_ox_cys = 0.050;
  const double k_drift = 0.0015;
  const double k_tet = 0.012;
  const double eta_d = 0.150;
  const double beta_b = 0.250;
  const double beta_h = 0.650;
  const double gamma_c = 0.340;
  const double theta_gsh = 1.150;
  const double theta_cys = 0.850;
  const double c1 = 0.220;
  const double c2 = 0.040;

  double W  = x[0] > 0.0 ? x[0] : 0.0;
  double D  = x[1] > 0.0 ? x[1] : 0.0;
  double B  = x[2] > 0.0 ? x[2] : 0.0;
  double E  = x[3] > 0.0 ? x[3] : 0.0;
  double G  = x[4] > 0.0 ? x[4] : 0.0;
  double N  = x[5] > 0.0 ? x[5] : 0.0;
  double C  = x[6] > 0.0 ? x[6] : 0.0;
  double M8 = x[7] > 0.0 ? x[7] : 0.0;

  double M = W + D;
  double h = (M > 1e-6) ? (D / M) : 0.0;

  double num_er = 1.0 + beta_b * B + beta_h * h;
  double eff_mass = W + eta_d * D;
  if (eff_mass < 0.10) eff_mass = 0.10;
  double den_er = eff_mass * (1.0 + gamma_c * C);
  double e_r = num_er / den_er;

  double f_eff = 1.0 - c1 * u1 + c2 * u2;
  double ep = e_r * f_eff * f_eff;

  double n_sat = N / (1.0 + N);
  double k_bio = k_bio0 + (k_biomax - k_bio0) * u2 * n_sat;
  double v_lys = K_B / (K_B + B);
  double k_mito_d = k_mito_w + (k_mito_d_max - k_mito_w) * u1;

  double cap_term = 1.0 - (M / K_cap);
  double ep_excess = (ep > theta_gsh) ? (ep - theta_gsh) : 0.0;
  double dmg_rate = k_dmg * ep_excess;

  const double b_basal = 472.0 / 65536.0;
  dxdt[0] = k_bio * W * cap_term - (k_mito_w + dmg_rate) * W;
  dxdt[1] = k_bio * k_rep_adv * D * cap_term - k_mito_d * v_lys * D + dmg_rate * W;
  dxdt[2] = b_basal + k_acc * D * (1.0 - u0) - k_clr * u0 * B;
  dxdt[3] = k_E * (ep - E);
  dxdt[4] = k_g_prod * E - k_g_clr * G;

  double e_excess = (E > theta_gsh) ? (E - theta_gsh) : 0.0;
  dxdt[5] = (k_n_syn0 + k_n_syn2 * u2) - k_parp * e_excess * N - k_n_use * N;

  double e_cys = (E > theta_cys) ? (E - theta_cys) : 0.0;
  dxdt[6] = (k_oct0 + k_oct2 * u2) * (1.30 - C) - k_ox_cys * e_cys * C;
  dxdt[7] = k_drift * E * (10.0 - M8) - k_tet * u2 * N * C * M8;
}

static void rk4_core8_step_fp64(double x[8], double u0, double u1, double u2, double dt) {
  double k1[8], k2[8], k3[8], k4[8], tmp[8];
  eval_core8_rhs_fp64(x, u0, u1, u2, k1);
  for (int i = 0; i < 8; ++i) tmp[i] = x[i] + 0.5 * dt * k1[i];
  eval_core8_rhs_fp64(tmp, u0, u1, u2, k2);
  for (int i = 0; i < 8; ++i) tmp[i] = x[i] + 0.5 * dt * k2[i];
  eval_core8_rhs_fp64(tmp, u0, u1, u2, k3);
  for (int i = 0; i < 8; ++i) tmp[i] = x[i] + dt * k3[i];
  eval_core8_rhs_fp64(tmp, u0, u1, u2, k4);
  for (int i = 0; i < 8; ++i) {
    x[i] += (dt / 6.0) * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]);
  }
}

/* Scalar RV32F IEEE-754 Single-Precision (`float`) 8-State ODE RHS & RK4 Stepper
 * (Evaluates scalar FPU accuracy on RV32IMF_Zve32x, where Zve32x lacks vector FP32) */
static void eval_core8_rhs_fp32(
    const float x[8],
    float u0,
    float u1,
    float u2,
    float dxdt[8]) {
  const float K_cap = 1.250f;
  const float k_rep_adv = 1.320f;
  const float k_bio0 = 0.075f;
  const float k_biomax = 0.140f;
  const float k_mito_w = 0.015f;
  const float k_mito_d_max = 0.115f;
  const float k_dmg = 0.040f;
  const float k_acc = 0.080f;
  const float k_clr = 0.180f;
  const float K_B = 0.250f;
  const float k_E = 0.250f;
  const float k_g_prod = 0.080f;
  const float k_g_clr = 0.160f;
  const float k_n_syn0 = 0.060f;
  const float k_n_syn2 = 0.080f;
  const float k_parp = 0.060f;
  const float k_n_use = 0.095f;
  const float k_oct0 = 0.040f;
  const float k_oct2 = 0.100f;
  const float k_ox_cys = 0.050f;
  const float k_drift = 0.0015f;
  const float k_tet = 0.012f;
  const float eta_d = 0.150f;
  const float beta_b = 0.250f;
  const float beta_h = 0.650f;
  const float gamma_c = 0.340f;
  const float theta_gsh = 1.150f;
  const float theta_cys = 0.850f;
  const float c1 = 0.220f;
  const float c2 = 0.040f;
  const float b_basal = 472.0f / 65536.0f;

  float W  = x[0] > 0.0f ? x[0] : 0.0f;
  float D  = x[1] > 0.0f ? x[1] : 0.0f;
  float B  = x[2] > 0.0f ? x[2] : 0.0f;
  float E  = x[3] > 0.0f ? x[3] : 0.0f;
  float G  = x[4] > 0.0f ? x[4] : 0.0f;
  float N  = x[5] > 0.0f ? x[5] : 0.0f;
  float C  = x[6] > 0.0f ? x[6] : 0.0f;
  float M8 = x[7] > 0.0f ? x[7] : 0.0f;

  float M = W + D;
  float h = (M > 1e-6f) ? (D / M) : 0.0f;

  float num_er = 1.0f + beta_b * B + beta_h * h;
  float eff_mass = W + eta_d * D;
  if (eff_mass < 0.10f) eff_mass = 0.10f;
  float den_er = eff_mass * (1.0f + gamma_c * C);
  float e_r = num_er / den_er;

  float f_eff = 1.0f - c1 * u1 + c2 * u2;
  float ep = e_r * f_eff * f_eff;

  float n_sat = N / (1.0f + N);
  float k_bio = k_bio0 + (k_biomax - k_bio0) * u2 * n_sat;
  float v_lys = K_B / (K_B + B);
  float k_mito_d = k_mito_w + (k_mito_d_max - k_mito_w) * u1;

  float cap_term = 1.0f - (M / K_cap);
  float ep_excess = (ep > theta_gsh) ? (ep - theta_gsh) : 0.0f;
  float dmg_rate = k_dmg * ep_excess;

  dxdt[0] = k_bio * W * cap_term - (k_mito_w + dmg_rate) * W;
  dxdt[1] = k_bio * k_rep_adv * D * cap_term - k_mito_d * v_lys * D + dmg_rate * W;
  dxdt[2] = b_basal + k_acc * D * (1.0f - u0) - k_clr * u0 * B;
  dxdt[3] = k_E * (ep - E);
  dxdt[4] = k_g_prod * E - k_g_clr * G;

  float e_excess = (E > theta_gsh) ? (E - theta_gsh) : 0.0f;
  dxdt[5] = (k_n_syn0 + k_n_syn2 * u2) - k_parp * e_excess * N - k_n_use * N;

  float e_cys = (E > theta_cys) ? (E - theta_cys) : 0.0f;
  dxdt[6] = (k_oct0 + k_oct2 * u2) * (1.30f - C) - k_ox_cys * e_cys * C;
  dxdt[7] = k_drift * E * (10.0f - M8) - k_tet * u2 * N * C * M8;
}

static void rk4_core8_step_fp32(float x[8], float u0, float u1, float u2, float dt) {
  float k1[8], k2[8], k3[8], k4[8], tmp[8];
  eval_core8_rhs_fp32(x, u0, u1, u2, k1);
  for (int i = 0; i < 8; ++i) tmp[i] = x[i] + 0.5f * dt * k1[i];
  eval_core8_rhs_fp32(tmp, u0, u1, u2, k2);
  for (int i = 0; i < 8; ++i) tmp[i] = x[i] + 0.5f * dt * k2[i];
  eval_core8_rhs_fp32(tmp, u0, u1, u2, k3);
  for (int i = 0; i < 8; ++i) tmp[i] = x[i] + dt * k3[i];
  eval_core8_rhs_fp32(tmp, u0, u1, u2, k4);
  for (int i = 0; i < 8; ++i) {
    x[i] += (dt / 6.0f) * (k1[i] + 2.0f * k2[i] + 2.0f * k3[i] + k4[i]);
  }
}

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

  /* 4. Test Rolling-Window Noise Filter + Secondary Modality Enforcement */
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
  /* Verify secondary-modality enforcement: if brood temp is 0 (outside [3330, 3570]), must NOT trigger! */
  mitorestorer_gate_action_t g_sec_fail = coralnpu_rolling_fuzzy_sensor_gate(
      CORAL_TARGET_APIS_SMART_HIVE, &guard_state, 2215, 0, 1080, true);
  assert(g_sec_fail == BCR_GATE_STAGE1_MITOPHAGY_PURGE);

  /* 5. Test CAR-T Amperometric Biofouling Fallback & Human Post-Meal Ketone Skew */
  mitorestorer_rolling_guard_dtcm_t cart_guard = {0};
  cart_guard.prev_lac_pyr_x100 = 2400;
  (void)coralnpu_rolling_fuzzy_sensor_gate(CORAL_TARGET_CART_GMP_BIOREACTOR, &cart_guard, 1600, 720, 420, false);
  (void)coralnpu_rolling_fuzzy_sensor_gate(CORAL_TARGET_CART_GMP_BIOREACTOR, &cart_guard, 2300, 720, 430, false);
  mitorestorer_gate_action_t g_biofoul_hold = coralnpu_rolling_fuzzy_sensor_gate(
      CORAL_TARGET_CART_GMP_BIOREACTOR, &cart_guard, 1100, 720, 440, false);
  assert(g_biofoul_hold == BCR_GATE_BIOFOUL_TIMED64H_OPTICAL_PH_HOLD);
  mitorestorer_gate_action_t g_biofoul_pass = coralnpu_rolling_fuzzy_sensor_gate(
      CORAL_TARGET_CART_GMP_BIOREACTOR, &cart_guard, 1100, 720, 640, true);
  assert(g_biofoul_pass == BCR_GATE_STAGE2_ANABOLIC_BUILD);

  mitorestorer_gate_action_t g_meal_skew = coralnpu_rolling_fuzzy_sensor_gate(
      CORAL_TARGET_HUMAN_WEARABLE_PATCH, &guard_state, 15, 2400, 360, true);
  assert(g_meal_skew == BCR_GATE_STAGE1_MITOPHAGY_PURGE);

  /* 6. Test Bidirectional Actuator Rollback & Reference Anchor Restoration */
  coralnpu_save_breadcrumb_ckpt(&guard_state, &hive_state_d4_5);
  int32_t saved_w_q16 = hive_state_d4_5.w_mass_q16;
  hive_state_d4_5.w_mass_q16 = (int32_t)(0.31 * 65536);

  mitorestorer_gate_action_t g_apis_rb = coralnpu_bidirectional_rollback_guard(
      CORAL_TARGET_APIS_SMART_HIVE, BCR_GATE_STAGE2_ANABOLIC_BUILD, &hive_state_d4_5, &guard_state, 3280, 72);
  assert(g_apis_rb == BCR_GATE_STAGE2_ROLLBACK_ABORT);
  assert(hive_state_d4_5.w_mass_q16 == saved_w_q16);

  mitorestorer_gate_action_t g_cart_rb = coralnpu_bidirectional_rollback_guard(
      CORAL_TARGET_CART_GMP_BIOREACTOR, BCR_GATE_STAGE2_ANABOLIC_BUILD, &cart_state_h64, &guard_state, 682, 4100);
  assert(g_cart_rb == BCR_GATE_STAGE2_ROLLBACK_ABORT);

  mitorestorer_gate_action_t g_apoe4_rb = coralnpu_bidirectional_rollback_guard(
      CORAL_TARGET_HUMAN_WEARABLE_PATCH, BCR_GATE_STAGE2_ANABOLIC_BUILD, &cart_state_h64, &guard_state, 4000, 680);
  assert(g_apoe4_rb == BCR_GATE_STAGE2_ROLLBACK_ABORT);

  /* 7. Full 1,000-Step 8-State 4-Stage RK4 Q16.16 vs. IEEE-754 FP64 Parity Verification */
  mitorestorer_core_params_q16_t p_q16 = {
    .k_cap_q16        = Q16_FROM_FLOAT(1.250),
    .k_rep_adv_q16    = Q16_FROM_FLOAT(1.320),
    .k_bio0_q16       = Q16_FROM_FLOAT(0.075),
    .k_biomax_q16     = Q16_FROM_FLOAT(0.140),
    .k_mito_w_q16     = Q16_FROM_FLOAT(0.015),
    .k_mito_d_max_q16 = Q16_FROM_FLOAT(0.115),
    .k_dmg_q16        = Q16_FROM_FLOAT(0.040),
    .k_acc_q16        = Q16_FROM_FLOAT(0.080),
    .k_clr_q16        = Q16_FROM_FLOAT(0.180),
    .k_b_half_q16     = Q16_FROM_FLOAT(0.250),
    .k_e_relax_q16    = Q16_FROM_FLOAT(0.250),
    .k_g_prod_q16     = Q16_FROM_FLOAT(0.080),
    .k_g_clr_q16      = Q16_FROM_FLOAT(0.160),
    .k_n_syn0_q16     = Q16_FROM_FLOAT(0.060),
    .k_n_syn2_q16     = Q16_FROM_FLOAT(0.080),
    .k_parp_q16       = Q16_FROM_FLOAT(0.060),
    .k_n_use_q16      = Q16_FROM_FLOAT(0.095),
    .k_oct0_q16       = Q16_FROM_FLOAT(0.040),
    .k_oct2_q16       = Q16_FROM_FLOAT(0.100),
    .k_ox_cys_q16     = Q16_FROM_FLOAT(0.050),
    .k_drift_q16      = Q16_FROM_FLOAT(0.0015),
    .k_tet_q16        = Q16_FROM_FLOAT(0.012),
    .eta_d_q16        = Q16_FROM_FLOAT(0.150),
    .beta_b_q16       = Q16_FROM_FLOAT(0.250),
    .beta_h_q16       = Q16_FROM_FLOAT(0.650),
    .gamma_c_q16      = Q16_FROM_FLOAT(0.340),
    .theta_gsh_q16    = Q16_FROM_FLOAT(1.150),
    .theta_cys_q16    = Q16_FROM_FLOAT(0.850),
    .c1_flux_q16      = Q16_FROM_FLOAT(0.220),
    .c2_flux_q16      = Q16_FROM_FLOAT(0.040),
    .h_gate_q16       = Q16_FROM_FLOAT(0.080),
    .b_gate_q16       = Q16_FROM_FLOAT(0.120)
  };

  double x_fp64[8] = {0.60, 0.40, 0.50, 1.35, 0.85, 0.70, 0.65, 4.50};
  mitorestorer_ode_state_q16_t x_q16 = {
    .w_mass_q16       = Q16_FROM_FLOAT(0.60),
    .d_mass_q16       = Q16_FROM_FLOAT(0.40),
    .bottleneck_b_q16 = Q16_FROM_FLOAT(0.50),
    .energy_pot_e_q16 = Q16_FROM_FLOAT(1.35),
    .stress_g_ng_q16  = Q16_FROM_FLOAT(0.85),
    .nad_pool_n_q16   = Q16_FROM_FLOAT(0.70),
    .pcr_buffer_c_q16 = Q16_FROM_FLOAT(0.65),
    .epigenetic_m8_q16= Q16_FROM_FLOAT(4.50)
  };

  int32_t dt_q16 = Q16_FROM_FLOAT(0.09); /* 5898 -> dt = 5898 / 65536 = 0.089996337890625 d */
  double dt_fp64 = (double)dt_q16 / 65536.0;
  double max_rel_err_pct[8] = {0};
  double sum_sq_rel_err[8] = {0};
  double global_max_err_pct = 0.0;

  bool stage2_latched = false;
  int switch_step = -1;
  int64_t rk4_rem[8] = {0};
  for (int step = 0; step < 1000; ++step) {
    int32_t* s_q16 = (int32_t*)&x_q16;
    if (!stage2_latched && coralnpu_core8_stage_gate_q16(&x_q16, p_q16.h_gate_q16, p_q16.b_gate_q16)) {
      stage2_latched = true;
      switch_step = step;
    }

    int32_t u0_q16 = Q16_ONE;
    int32_t u1_q16 = stage2_latched ? Q16_FROM_FLOAT(0.35) : Q16_ONE; /* 0.35 maintenance in Stage 2 */
    int32_t u2_q16 = stage2_latched ? Q16_ONE : 0;

    double u0_f = 1.0;
    double u1_f = stage2_latched ? (double)Q16_FROM_FLOAT(0.35) / 65536.0 : 1.0;
    double u2_f = stage2_latched ? 1.0 : 0.0;

    coralnpu_rk4_core8_full_step_q16(&x_q16, &p_q16, u0_q16, u1_q16, u2_q16, dt_q16, rk4_rem);
    rk4_core8_step_fp64(x_fp64, u0_f, u1_f, u2_f, dt_fp64);

    for (int i = 0; i < 8; ++i) {
      double q_val = Q16_TO_DOUBLE(s_q16[i]);
      double ref = fabs(x_fp64[i]) > 1e-6 ? fabs(x_fp64[i]) : 1e-6;
      double rel_pct = 100.0 * fabs(q_val - x_fp64[i]) / ref;
      if (rel_pct > max_rel_err_pct[i]) max_rel_err_pct[i] = rel_pct;
      if (rel_pct > global_max_err_pct) global_max_err_pct = rel_pct;
      sum_sq_rel_err[i] += rel_pct * rel_pct;
    }
  }
  printf("CORAL_NPU_C99_HIL_PASS: Apis_d0=%d Apis_d4_5=%d CART_h64=%d Human_Override=%d WindGust_Rej=%d Biofoul_Hold=%d MealSkew_Rej=%d Rollbacks=%d/%d/%d Logit=%d SkippedBins=%d/16\n",
         g_apis_d0, g_apis_d4_5, g_cart_h64, g_hum_override, g_wind_gust, g_biofoul_hold, g_meal_skew,
         g_apis_rb, g_cart_rb, g_apoe4_rb, logit_acc, skipped_bins);

  int32_t* final_q16 = (int32_t*)&x_q16;
  for (int i = 0; i < 8; ++i) {
    double rms_pct = sqrt(sum_sq_rel_err[i] / 1000.0);
    printf("RK4_PARITY_STATE[%d]: FP64=%.6f | Q16_INT=%d (%.6f) | TermErr=%.4f%% | MaxErr=%.4f%% | RMSErr=%.4f%%\n",
           i, x_fp64[i], final_q16[i], Q16_TO_DOUBLE(final_q16[i]),
           100.0 * fabs(Q16_TO_DOUBLE(final_q16[i]) - x_fp64[i]) / (fabs(x_fp64[i]) > 1e-3 ? fabs(x_fp64[i]) : 1e-3),
           max_rel_err_pct[i], rms_pct);
  }
  assert(global_max_err_pct < 0.50);
  (void)switch_step;

  /* 8. Dual-Path RVV Zve32x Widening Vector (`vwmul.vx` + `vnsra.wx ..., 16`) 0-ULP Parity Audit */
  int max_rvv_ulp_diff = 0;
  for (int trial = 0; trial < 256; ++trial) {
    mitorestorer_ode_state_q16_t st_scalar = {0};
    int32_t st_rvv_isa[8] = {0};
    int32_t k1_v[8], k2_v[8];
    for (int lane = 0; lane < 8; ++lane) {
      int32_t base_val = (int32_t)((trial * 131 + lane * 977 - 16000) * 19);
      ((int32_t*)&st_scalar)[lane] = base_val;
      st_rvv_isa[lane] = base_val;
      k1_v[lane] = (int32_t)((lane % 2 == 0 ? 1 : -1) * (trial * 37 + lane * 211));
      k2_v[lane] = (int32_t)((lane % 2 == 0 ? 1 : -1) * (trial * 41 + lane * 163));
    }
    coralnpu_rvv_rk4_bioenergetic_step(&st_scalar, k1_v, k2_v, dt_q16);
    coralnpu_zve32x_isa_emulated_step_q16(st_rvv_isa, k1_v, k2_v, dt_q16);
    for (int lane = 0; lane < 8; ++lane) {
      int diff = ((int32_t*)&st_scalar)[lane] - st_rvv_isa[lane];
      if (diff < 0) diff = -diff;
      if (diff > max_rvv_ulp_diff) max_rvv_ulp_diff = diff;
    }
  }
  assert(max_rvv_ulp_diff == 0);
  printf("RVV_ZVE32X_WIDENING_PARITY: Trials=256 | Lanes=8 | MaxULPDiff=%d\n", max_rvv_ulp_diff);

  /* 9. Multi-Step-Size Arithmetic & Fixed-Point Deadband Ablation (T_max = 90 d):
   *    Compares Scalar RV32F FP32 vs. Naive Q16.16 RK4 vs. Q8.24 No-Rem vs. Q8.24 + Sigma-Delta Rem */
  const int step_counts[5] = {500, 1000, 2000, 4000, 8000};
  for (int sc_idx = 0; sc_idx < 5; ++sc_idx) {
    int n_steps = step_counts[sc_idx];
    double dt_nom = 90.0 / (double)n_steps;
    int32_t dt_ab_q16 = Q16_FROM_FLOAT(dt_nom);
    double dt_ab_fp64 = (double)dt_ab_q16 / 65536.0;
    float dt_ab_fp32 = (float)dt_ab_fp64;

    double x_ref[8] = {0.60, 0.40, 0.50, 1.35, 0.85, 0.70, 0.65, 4.50};
    float  x_f32[8] = {0.60f, 0.40f, 0.50f, 1.35f, 0.85f, 0.70f, 0.65f, 4.50f};
    mitorestorer_ode_state_q16_t x_naive = {
      .w_mass_q16 = Q16_FROM_FLOAT(0.60), .d_mass_q16 = Q16_FROM_FLOAT(0.40),
      .bottleneck_b_q16 = Q16_FROM_FLOAT(0.50), .energy_pot_e_q16 = Q16_FROM_FLOAT(1.35),
      .stress_g_ng_q16 = Q16_FROM_FLOAT(0.85), .nad_pool_n_q16 = Q16_FROM_FLOAT(0.70),
      .pcr_buffer_c_q16 = Q16_FROM_FLOAT(0.65), .epigenetic_m8_q16 = Q16_FROM_FLOAT(4.50)
    };
    mitorestorer_ode_state_q16_t x_norem = x_naive;
    mitorestorer_ode_state_q16_t x_rem   = x_naive;
    int64_t rem_acc[8] = {0};

    bool s2_ref = false, s2_f32 = false, s2_naive = false, s2_norem = false, s2_rem = false;
    double max_err_f32 = 0.0, max_err_naive = 0.0, max_err_norem = 0.0, max_err_rem = 0.0;

    for (int st = 0; st < n_steps; ++st) {
      if (!s2_ref && (x_ref[1] / (x_ref[0] + x_ref[1])) <= 0.080 && x_ref[2] <= 0.120) s2_ref = true;
      if (!s2_f32 && (x_f32[1] / (x_f32[0] + x_f32[1])) <= 0.080f && x_f32[2] <= 0.120f) s2_f32 = true;
      if (!s2_naive && coralnpu_core8_stage_gate_q16(&x_naive, p_q16.h_gate_q16, p_q16.b_gate_q16)) s2_naive = true;
      if (!s2_norem && coralnpu_core8_stage_gate_q16(&x_norem, p_q16.h_gate_q16, p_q16.b_gate_q16)) s2_norem = true;
      if (!s2_rem   && coralnpu_core8_stage_gate_q16(&x_rem,   p_q16.h_gate_q16, p_q16.b_gate_q16)) s2_rem   = true;

      double u1_r = s2_ref ? (double)Q16_FROM_FLOAT(0.35) / 65536.0 : 1.0;
      double u2_r = s2_ref ? 1.0 : 0.0;
      rk4_core8_step_fp64(x_ref, 1.0, u1_r, u2_r, dt_ab_fp64);

      float u1_32 = s2_f32 ? (float)Q16_FROM_FLOAT(0.35) / 65536.0f : 1.0f;
      float u2_32 = s2_f32 ? 1.0f : 0.0f;
      rk4_core8_step_fp32(x_f32, 1.0f, u1_32, u2_32, dt_ab_fp32);

      coralnpu_rk4_core8_naive_q16_step(
          &x_naive, &p_q16, Q16_ONE,
          s2_naive ? Q16_FROM_FLOAT(0.35) : Q16_ONE,
          s2_naive ? Q16_ONE : 0, dt_ab_q16);
      coralnpu_rk4_core8_full_step_q16(
          &x_norem, &p_q16, Q16_ONE,
          s2_norem ? Q16_FROM_FLOAT(0.35) : Q16_ONE,
          s2_norem ? Q16_ONE : 0, dt_ab_q16, NULL);
      coralnpu_rk4_core8_full_step_q16(
          &x_rem, &p_q16, Q16_ONE,
          s2_rem ? Q16_FROM_FLOAT(0.35) : Q16_ONE,
          s2_rem ? Q16_ONE : 0, dt_ab_q16, rem_acc);

      for (int i = 0; i < 8; ++i) {
        double denom_r = fabs(x_ref[i]) > 1e-4 ? fabs(x_ref[i]) : 1e-4;
        double e_f32   = 100.0 * fabs((double)x_f32[i] - x_ref[i]) / denom_r;
        double e_naive = 100.0 * fabs(Q16_TO_DOUBLE(((int32_t*)&x_naive)[i]) - x_ref[i]) / denom_r;
        double e_norem = 100.0 * fabs(Q16_TO_DOUBLE(((int32_t*)&x_norem)[i]) - x_ref[i]) / denom_r;
        double e_rem   = 100.0 * fabs(Q16_TO_DOUBLE(((int32_t*)&x_rem)[i])   - x_ref[i]) / denom_r;
        if (e_f32   > max_err_f32)   max_err_f32   = e_f32;
        if (e_naive > max_err_naive) max_err_naive = e_naive;
        if (e_norem > max_err_norem) max_err_norem = e_norem;
        if (e_rem   > max_err_rem)   max_err_rem   = e_rem;
      }
    }
    printf("FORMAT_ABLATION[N=%d,dt=%.5f]: FP32_MaxErr=%.4f%% | NaiveQ16_MaxErr=%.4f%% | Q24NoRem_MaxErr=%.4f%% | Q24Rem_MaxErr=%.4f%%\n",
           n_steps, dt_ab_fp64, max_err_f32, max_err_naive, max_err_norem, max_err_rem);
  }

  /* 10. Bare-Metal C99 Q16.16 7-State Huber-EKF Closed-Loop Test & Compiler DTCM sizeof() Audit */
  coralnpu_dtcm_bank0_workspace_t dtcm_ws = {0};
  size_t dtcm_bytes = sizeof(coralnpu_dtcm_bank0_workspace_t);
  assert(dtcm_bytes <= 8192);

  double x_plant[8] = {0.58, 0.42, 0.52, 1.38, 0.88, 0.68, 0.63, 4.50}; /* Perturbed true plant */
  dtcm_ws.params_q16 = p_q16;
  dtcm_ws.x_hat_q16.w_mass_q16        = Q16_FROM_FLOAT(0.60);
  dtcm_ws.x_hat_q16.d_mass_q16        = Q16_FROM_FLOAT(0.40);
  dtcm_ws.x_hat_q16.bottleneck_b_q16  = Q16_FROM_FLOAT(0.50);
  dtcm_ws.x_hat_q16.energy_pot_e_q16  = Q16_FROM_FLOAT(1.35);
  dtcm_ws.x_hat_q16.stress_g_ng_q16   = Q16_FROM_FLOAT(0.85);
  dtcm_ws.x_hat_q16.nad_pool_n_q16    = Q16_FROM_FLOAT(0.70);
  dtcm_ws.x_hat_q16.pcr_buffer_c_q16  = Q16_FROM_FLOAT(0.65);
  dtcm_ws.x_hat_q16.epigenetic_m8_q16 = Q16_FROM_FLOAT(4.50);

  const double p_init_diag[7] = {0.02, 0.02, 0.03, 0.04, 0.03, 0.02, 0.02};
  const double q_rate_diag[7] = {1e-4, 1e-4, 2e-4, 3e-4, 2e-4, 1e-4, 1e-4};
  for (int i = 0; i < 7; ++i) {
    dtcm_ws.ekf_P_cov_q24[i][i] = (int32_t)(p_init_diag[i] * 16777216.0 + 0.5);
    dtcm_ws.ekf_Q_diag_q24[i]   = (int32_t)(q_rate_diag[i] * dt_fp64 * 16777216.0 + 0.5);
    if (dtcm_ws.ekf_Q_diag_q24[i] < 1) dtcm_ws.ekf_Q_diag_q24[i] = 1;
  }
  int32_t r_var_q24 = (int32_t)(0.030 * 0.030 * 16777216.0 + 0.5); /* sigma_v = 0.03 -> R = 0.0009 (15,099 in Q8.24) */
  int32_t gamma_huber_q16 = Q16_FROM_FLOAT(1.50);
  int32_t h_gate_margin_q16 = Q16_FROM_FLOAT(0.065);

  uint32_t lcg = 20261109u;
  double ar_v[4] = {0.0, 0.0, 0.0, 0.0};
  bool ekf_s2 = false;
  int ekf_persist = 0;
  int ekf_sw_step = -1;
  double sum_sq_h_err = 0.0;

  for (int st = 0; st < 1000; ++st) {
    /* Deterministic Box-Muller / Irwin-Hall Gaussian AR(1) (rho=0.60) + 5% 6.5-sigma spikes */
    int32_t y_meas_q16[4];
    double atp_true = (x_plant[0] + 0.150 * x_plant[1]) * (1.0 + 0.340 * x_plant[6]);
    double y_clean[4] = {x_plant[3], atp_true, x_plant[4], x_plant[5]};
    for (int ch = 0; ch < 4; ++ch) {
      double u_sum = 0.0;
      for (int m = 0; m < 12; ++m) {
        lcg = lcg * 1664525u + 1013904223u;
        u_sum += (double)(lcg >> 8) / 16777216.0;
      }
      double z_norm = u_sum - 6.0;
      ar_v[ch] = 0.60 * ar_v[ch] + 0.80 * 0.030 * z_norm;
      lcg = lcg * 1664525u + 1013904223u;
      double u_sp = (double)(lcg >> 8) / 16777216.0;
      double spike = (u_sp < 0.025) ? (6.5 * 0.030) : ((u_sp > 0.975) ? (-6.5 * 0.030) : 0.0);
      y_meas_q16[ch] = Q16_FROM_FLOAT(y_clean[ch] + ar_v[ch] + spike);
    }

    if (st > 0) {
      coralnpu_huber_ekf7_meas_update_q16(
          &dtcm_ws.x_hat_q16, dtcm_ws.ekf_P_cov_q24, &dtcm_ws.params_q16,
          y_meas_q16, r_var_q24, gamma_huber_q16);
    }

    if (!ekf_s2) {
      if (coralnpu_core8_stage_gate_q16(&dtcm_ws.x_hat_q16, h_gate_margin_q16, dtcm_ws.params_q16.b_gate_q16)) {
        ekf_persist++;
        if (ekf_persist >= 8) {
          ekf_s2 = true;
          ekf_sw_step = st;
        }
      } else {
        ekf_persist = 0;
      }
    }

    double h_tr = x_plant[1] / (x_plant[0] + x_plant[1]);
    double h_eh = Q16_TO_DOUBLE(dtcm_ws.x_hat_q16.d_mass_q16) /
                  (Q16_TO_DOUBLE(dtcm_ws.x_hat_q16.w_mass_q16) + Q16_TO_DOUBLE(dtcm_ws.x_hat_q16.d_mass_q16));
    sum_sq_h_err += (h_eh - h_tr) * (h_eh - h_tr);

    int32_t u1_c = ekf_s2 ? Q16_FROM_FLOAT(0.35) : Q16_ONE;
    int32_t u2_c = ekf_s2 ? Q16_ONE : 0;
    rk4_core8_step_fp64(x_plant, 1.0, ekf_s2 ? (double)Q16_FROM_FLOAT(0.35) / 65536.0 : 1.0, ekf_s2 ? 1.0 : 0.0, dt_fp64);
    coralnpu_huber_ekf7_predict_q16(
        &dtcm_ws.x_hat_q16, dtcm_ws.ekf_P_cov_q24, &dtcm_ws.params_q16,
        Q16_ONE, u1_c, u2_c, dt_q16, dtcm_ws.ekf_Q_diag_q24, dtcm_ws.rk4_sub_lsb_rem);
  }
  double ekf_h_rmse_pct = 100.0 * sqrt(sum_sq_h_err / 1000.0);
  double ekf_h_term = x_plant[1] / (x_plant[0] + x_plant[1]);
  assert(ekf_s2 == true && ekf_h_rmse_pct < 3.0 && x_plant[0] > 0.95);
  printf("C99_HUBER_EKF_CLOSED_LOOP: Steps=1000 | SwitchStep=%d | h_RMSE=%.2f%% | W_term=%.4f | h_term=%.4f | DTCM_Workspace=%zuB/8192B\n",
         ekf_sw_step, ekf_h_rmse_pct, x_plant[0], ekf_h_term, dtcm_bytes);

  return 0;
}

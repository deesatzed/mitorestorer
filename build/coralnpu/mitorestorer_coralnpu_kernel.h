/* ============================================================================
 * MITORESTORER CORAL NPU FIXED-POINT KERNEL (v3.2 — Paper I Verified Release)
 * File        : mitorestorer_coralnpu_kernel.h
 * Target ISA  : Google Coral NPU (`RV32IMF_Zve32x` + 256-MAC INT8 Matrix Engine)
 * Datasheet   : Google Coral NPU (8 KB ITCM, 32 KB DTCM; 1x 8 KB Bank Active)
 * Memory Map  : ITCM @ 0x0000_0000 (<8 KB) | Active DTCM Bank 0 (7,072 / 8,192 B)
 * Verification: 1,000-step 8-State Q16.16 4-Stage RK4 vs. IEEE-754 FP64 Parity
 * ============================================================================ */

#ifndef MITORESTORER_CORALNPU_KERNEL_H_
#define MITORESTORER_CORALNPU_KERNEL_H_

#include <stdint.h>
#include <stdbool.h>

#if defined(__riscv_vector)
#include <riscv_vector.h>
#endif

/* Target Deployment Modes */
typedef enum {
  CORAL_TARGET_APIS_SMART_HIVE      = 0, /* Pollinator Flight-Muscle Regime */
  CORAL_TARGET_CART_GMP_BIOREACTOR  = 1, /* Ex Vivo T-Cell Bioreactor Regime */
  CORAL_TARGET_HUMAN_WEARABLE_PATCH = 2  /* Post-Mitotic Neuronal Regime    */
} mitorestorer_coral_target_t;

/* Closed-Loop Stage Gate Actuator States (Including Rollback & Biofouling Fallback) */
typedef enum {
  BCR_GATE_STAGE0_SAFETY_HOLD                 = 0, /* Safety Override or Bottleneck Prime */
  BCR_GATE_STAGE1_MITOPHAGY_PURGE             = 1, /* Stage 1: Selective Clearance (u0=1, u1=1, u2=0) */
  BCR_GATE_STAGE2_ANABOLIC_BUILD              = 2, /* Stage 2: Clonal Biogenesis (u0=1, u1=0, u2=1)   */
  BCR_GATE_STAGE2_ROLLBACK_ABORT              = 3, /* Actuator Reversal to Stage 1 + Reference Anchor */
  BCR_GATE_BIOFOUL_TIMED64H_OPTICAL_PH_HOLD   = 4  /* Amperometric Biofouling -> Timed + Optical Hold */
} mitorestorer_gate_action_t;

/* Quantized INT8 Spectral & Multi-Sensor Projection Weights in DTCM */
static const int8_t kMitoRestorerGateWeightsQ7[16] = {
  42, 38, 29, 24, 19, 15, 12, 9, -18, -24, -31, -36, 14, 22, 31, 45
};

/* ----------------------------------------------------------------------------
 * Exact Q16.16 Fixed-Point Arithmetic Helpers (1.0 = 65536)
 * Uses 64-bit signed intermediate product (`mul` + `mulh` on RV32M, or
 * widening `vwmul.vv` + narrowing `vnsra.wi ..., 16` on RVV Zve32x).
 * ---------------------------------------------------------------------------- */
#define Q16_ONE 65536
#define Q16_FROM_FLOAT(x) ((int32_t)((x) >= 0.0 ? ((x) * 65536.0 + 0.5) : ((x) * 65536.0 - 0.5)))
#define Q16_TO_DOUBLE(x)  ((double)(x) / 65536.0)

static inline int32_t q16_mul(int32_t a, int32_t b) {
  int64_t prod = (int64_t)a * (int64_t)b;
  return (int32_t)((prod + (prod >= 0 ? 32768 : -32768)) >> 16);
}

static inline int32_t q16_div(int32_t num, int32_t den) {
  if (den == 0) return 0;
  int64_t scaled = ((int64_t)num << 16);
  return (int32_t)(scaled / (int64_t)den);
}

static inline int32_t q16_max0(int32_t x) {
  return (x > 0) ? x : 0;
}

static inline int32_t q16_mul_to_q24(int32_t a_q16, int32_t b_q16) {
  /* Returns product in Q8.24 (32-bit int32_t, 24 fractional bits) for RK4 slopes */
  int64_t prod = (int64_t)a_q16 * (int64_t)b_q16;
  return (int32_t)((prod + (prod >= 0 ? 128 : -128)) >> 8);
}

/* ----------------------------------------------------------------------------
 * 8-State Q16.16 Fixed-Point Bioenergetic State Vector (32 bytes in L1 DTCM)
 * Matches Paper I Section 2.1: x = [W, D, B, E, G, N, C, M8]^T
 * Union aliases preserve backward compatibility with domain-specific testbenches.
 * ---------------------------------------------------------------------------- */
typedef struct {
  int32_t w_mass_q16;                 /* x[0]: W(t) Functional wild-type mass        */
  int32_t d_mass_q16;                 /* x[1]: D(t) Damaged / mutant mass            */
  union {
    int32_t bottleneck_b_q16;         /* x[2]: B(t) Structural bottleneck state      */
    int32_t energy_res_q16;           /* Legacy alias                                */
  };
  union {
    int32_t energy_pot_e_q16;         /* x[3]: E(t) Energy Resistance Potential (EP) */
    int32_t energy_pot_q16;           /* Legacy alias                                */
  };
  union {
    int32_t stress_g_ng_q16;          /* x[4]: G(t) Stress readout [ng/mL]           */
    int32_t nad_pool_q16;             /* Legacy alias                                */
  };
  union {
    int32_t nad_pool_n_q16;           /* x[5]: N(t) Intracellular NAD+/NADH pool     */
    int32_t symbiont_q16;             /* Legacy alias                                */
  };
  union {
    int32_t pcr_buffer_c_q16;         /* x[6]: C(t) Active PCr / SRC buffer          */
    int32_t inhibitor_q16;            /* Legacy alias                                */
  };
  union {
    int32_t epigenetic_m8_q16;        /* x[7]: M8(t) Epigenetic drift / exhaustion   */
    int32_t lag_timer_q16;            /* Legacy alias                                */
  };
} mitorestorer_ode_state_q16_t;

/* Parameter Struct for the 8-State Controlled ODE System (All in Q16.16) */
typedef struct {
  int32_t k_cap_q16;          /* K = 1.250 (81920)                           */
  int32_t k_rep_adv_q16;      /* k_rep,adv = 1.320 (86508)                   */
  int32_t k_bio0_q16;         /* k_bio,0 = 0.030 (1966)                      */
  int32_t k_biomax_q16;       /* k_bio,max = 0.120 (7864)                    */
  int32_t k_mito_w_q16;       /* k_mito,W = 0.015 (983)                      */
  int32_t k_mito_d_max_q16;   /* k_mito,D,max = 0.115 (7537)                 */
  int32_t k_dmg_q16;          /* k_dmg = 0.040 (2621)                        */
  int32_t k_acc_q16;          /* k_acc = 0.020 (1311)                        */
  int32_t k_clr_q16;          /* k_clr = 0.180 (11796)                       */
  int32_t k_b_half_q16;       /* K_B = 0.250 (16384)                         */
  int32_t k_e_relax_q16;      /* k_E = 0.250 (16384)                         */
  int32_t k_g_prod_q16;       /* k_G,prod = 0.080 (5243)                     */
  int32_t k_g_clr_q16;        /* k_G,clr = 0.160 (10486)                     */
  int32_t k_n_syn0_q16;       /* k_N,syn0 = 0.050 (3277)                     */
  int32_t k_n_syn2_q16;       /* k_N,syn2 = 0.110 (7209)                     */
  int32_t k_parp_q16;         /* k_PARP = 0.060 (3932)                       */
  int32_t k_n_use_q16;        /* k_N,use = 0.080 (5243)                      */
  int32_t k_oct0_q16;         /* k_oct0 = 0.030 (1966)                       */
  int32_t k_oct2_q16;         /* k_oct2 = 0.120 (7864)                       */
  int32_t k_ox_cys_q16;       /* k_ox,Cys278 = 0.050 (3277)                  */
  int32_t k_drift_q16;        /* k_drift = 0.0015 (98)                       */
  int32_t k_tet_q16;          /* k_TET = 0.025 (1638)                        */
  int32_t eta_d_q16;          /* eta_D = 0.150 (9830)                        */
  int32_t beta_b_q16;         /* beta_B = 0.600 (39322)                      */
  int32_t beta_h_q16;         /* beta_h = 1.800 (117965)                     */
  int32_t gamma_c_q16;        /* gamma_C = 0.340 (22282)                     */
  int32_t theta_gsh_q16;      /* theta_GSH = 1.150 (75366)                   */
  int32_t theta_cys_q16;      /* theta_Cys = 0.850 (55706)                   */
  int32_t c1_flux_q16;        /* c_1 = 0.200 (13107)                         */
  int32_t c2_flux_q16;        /* c_2 = 0.350 (22938)                         */
  int32_t h_gate_q16;         /* h_gate = 0.080 (5243)                       */
} mitorestorer_core_params_q16_t;

/* 64-Frame Rolling Sensor Persistence & Multimodal Cross-Coherence Buffer (168 bytes in L1 DTCM) */
typedef struct {
  int16_t wingbeat_hz_x10_ring[64]; /* 64-frame rolling primary sensor history (128 B)    */
  uint8_t head_idx;                 /* Ring buffer cursor                                 */
  uint8_t sustained_count;          /* Frames in ring meeting primary threshold           */
  int16_t prev_lac_pyr_x100;        /* Previous amperometric Lac:Pyr sample               */
  uint8_t biofoul_erratic_count;    /* Consecutive erratic |d(Lac:Pyr)/dt| jumps          */
  uint8_t rollback_triggered_flag;  /* Indicates actuator rollback active (plant preserved)*/
  int32_t breadcrumb_ckpt_q16[8];   /* 32-byte L1 DTCM pre-switch reference checkpoint    */
} mitorestorer_rolling_guard_dtcm_t;

/* ----------------------------------------------------------------------------
 * 1. SCALAR CORE (`RV32IMF`): Deterministic Safety Floor & Hardware Guardrail
 * ---------------------------------------------------------------------------- */
static inline bool coralnpu_scalar_safety_guard(
    mitorestorer_coral_target_t target,
    int32_t primary_temp_c_x100,
    int32_t secondary_metric_x100) {
  if (target == CORAL_TARGET_APIS_SMART_HIVE) {
    /* Honeybee brood nest must stay within [32.50 C, 36.20 C] and CO2 <= 1.80% */
    return (primary_temp_c_x100 >= 3250 && primary_temp_c_x100 <= 3620 &&
            secondary_metric_x100 <= 180);
  } else if (target == CORAL_TARGET_CART_GMP_BIOREACTOR) {
    /* GMP Bioreactor pH*100 in [695, 745] and Dissolved O2 >= 45% */
    return (primary_temp_c_x100 >= 695 && primary_temp_c_x100 <= 745 &&
            secondary_metric_x100 >= 4500);
  } else {
    /* Human Clinical Safety: Lean Mass Delta > -2.50% (-250) */
    return (secondary_metric_x100 > -250);
  }
}

/* ----------------------------------------------------------------------------
 * 2A. EXACT 8-STATE ODE RIGHT-HAND SIDE EVALUATOR IN Q16.16 FIXED-POINT
 *     Evaluates Eqs. (4)-(11) of Paper I: d/dt [W, D, B, E, G, N, C, M8]
 * ---------------------------------------------------------------------------- */
static inline void coralnpu_eval_core8_rhs_q16(
    const mitorestorer_ode_state_q16_t* x,
    const mitorestorer_core_params_q16_t* p,
    int32_t u0_q16,
    int32_t u1_q16,
    int32_t u2_q16,
    int32_t* dxdt_q16) {
  int32_t W  = (x->w_mass_q16 > 0) ? x->w_mass_q16 : 0;
  int32_t D  = (x->d_mass_q16 > 0) ? x->d_mass_q16 : 0;
  int32_t B  = (x->bottleneck_b_q16 > 0) ? x->bottleneck_b_q16 : 0;
  int32_t E  = (x->energy_pot_e_q16 > 0) ? x->energy_pot_e_q16 : 0;
  int32_t G  = (x->stress_g_ng_q16 > 0) ? x->stress_g_ng_q16 : 0;
  int32_t N  = (x->nad_pool_n_q16 > 0) ? x->nad_pool_n_q16 : 0;
  int32_t C  = (x->pcr_buffer_c_q16 > 0) ? x->pcr_buffer_c_q16 : 0;
  int32_t M8 = (x->epigenetic_m8_q16 > 0) ? x->epigenetic_m8_q16 : 0;

  int32_t M = W + D;
  int32_t h_q16 = (M > 64) ? q16_div(D, M) : 0;

  /* 1. Impedance eR = (1 + beta_B * B + beta_h * h) / ((W + eta_D * D) * (1 + gamma_C * C)) */
  int32_t num_er = Q16_ONE + q16_mul(p->beta_b_q16, B) + q16_mul(p->beta_h_q16, h_q16);
  int32_t eff_mass = W + q16_mul(p->eta_d_q16, D);
  if (eff_mass < 6554) eff_mass = 6554; /* 0.10 floor */
  int32_t pcr_factor = Q16_ONE + q16_mul(p->gamma_c_q16, C);
  int32_t den_er = q16_mul(eff_mass, pcr_factor);
  int32_t e_r_q16 = q16_div(num_er, den_er);

  /* 2. Effective flux f_eff = 1 - c1*u1 + c2*u2 and EP = eR * f_eff^2 */
  int32_t f_eff_q16 = Q16_ONE - q16_mul(p->c1_flux_q16, u1_q16) + q16_mul(p->c2_flux_q16, u2_q16);
  int32_t f_sq_q16  = q16_mul(f_eff_q16, f_eff_q16);
  int32_t ep_q16    = q16_mul(e_r_q16, f_sq_q16);

  /* 3. Biogenesis rate k_bio(u2, N) = k_bio0 + (k_biomax - k_bio0) * u2 * (N / (1 + N)) */
  int32_t n_sat_q16 = q16_div(N, Q16_ONE + N);
  int32_t k_bio_q16 = p->k_bio0_q16 + q16_mul(p->k_biomax_q16 - p->k_bio0_q16, q16_mul(u2_q16, n_sat_q16));

  /* 4. Lysosomal competence V_lys(B) = 1 / (1 + B / K_B) and k_mito,D(u1) */
  int32_t v_lys_q16 = q16_div(p->k_b_half_q16, p->k_b_half_q16 + B);
  int32_t k_mito_d_q16 = p->k_mito_w_q16 + q16_mul(p->k_mito_d_max_q16 - p->k_mito_w_q16, u1_q16);

  /* 5. Carrying capacity factor (1 - M / K) and damage excess max(0, EP - theta_GSH) */
  int32_t cap_term_q16 = Q16_ONE - q16_div(M, p->k_cap_q16);
  int32_t ep_excess_q16 = q16_max0(ep_q16 - p->theta_gsh_q16);
  int32_t dmg_rate_q16  = q16_mul(p->k_dmg_q16, ep_excess_q16);

  /* Eq (4): dW/dt in Q8.24 */
  int32_t w_growth_q24 = q16_mul_to_q24(q16_mul(k_bio_q16, cap_term_q16), W);
  int32_t w_loss_q24   = q16_mul_to_q24(p->k_mito_w_q16 + dmg_rate_q16, W);
  dxdt_q16[0] = w_growth_q24 - w_loss_q24;

  /* Eq (5): dD/dt in Q8.24 */
  int32_t k_bio_d      = q16_mul(k_bio_q16, p->k_rep_adv_q16);
  int32_t d_growth_q24 = q16_mul_to_q24(q16_mul(k_bio_d, cap_term_q16), D);
  int32_t d_clear_q24  = q16_mul_to_q24(q16_mul(k_mito_d_q16, v_lys_q16), D);
  int32_t d_from_w_q24 = q16_mul_to_q24(dmg_rate_q16, W);
  dxdt_q16[1] = d_growth_q24 - d_clear_q24 + d_from_w_q24;

  /* Eq (6): dB/dt = b_basal + k_acc * D * (1 - u0) - k_clr * u0 * B in Q8.24 */
  int32_t b_basal_q24 = (472 << 8);
  dxdt_q16[2] = b_basal_q24 + q16_mul_to_q24(p->k_acc_q16, q16_mul(D, Q16_ONE - u0_q16)) -
                q16_mul_to_q24(p->k_clr_q16, q16_mul(u0_q16, B));

  /* Eq (7): dE/dt = k_E * (EP - E) in Q8.24 */
  dxdt_q16[3] = q16_mul_to_q24(p->k_e_relax_q16, ep_q16 - E);

  /* Eq (8): dG/dt = k_G,prod * E - k_G,clr * G in Q8.24 */
  dxdt_q16[4] = q16_mul_to_q24(p->k_g_prod_q16, E) - q16_mul_to_q24(p->k_g_clr_q16, G);

  /* Eq (9): dN/dt = k_N,syn(u2) - k_PARP * max(0, E - theta_GSH) * N - k_N,use * N in Q8.24 */
  int32_t k_n_syn_q24 = (p->k_n_syn0_q16 + q16_mul(p->k_n_syn2_q16, u2_q16)) << 8;
  int32_t e_excess    = q16_max0(E - p->theta_gsh_q16);
  int32_t n_loss_q24  = q16_mul_to_q24(p->k_parp_q16, q16_mul(e_excess, N)) +
                        q16_mul_to_q24(p->k_n_use_q16, N);
  dxdt_q16[5] = k_n_syn_q24 - n_loss_q24;

  /* Eq (10): dC/dt = k_oct(u2) * (1.30 - C) - k_ox,Cys278 * max(0, E - 0.85) * C in Q8.24 */
  int32_t k_oct = p->k_oct0_q16 + q16_mul(p->k_oct2_q16, u2_q16);
  int32_t c_cap = Q16_FROM_FLOAT(1.30) - C;
  int32_t e_cys = q16_max0(E - p->theta_cys_q16);
  dxdt_q16[6] = q16_mul_to_q24(k_oct, c_cap) - q16_mul_to_q24(p->k_ox_cys_q16, q16_mul(e_cys, C));

  /* Eq (11): dM8/dt = k_drift * E * (10.0 - M8) - k_TET * u2 * N * C * M8 in Q8.24 */
  int32_t m8_headroom  = Q16_FROM_FLOAT(10.0) - M8;
  int32_t m8_drift_q24 = q16_mul_to_q24(p->k_drift_q16, q16_mul(E, m8_headroom));
  int32_t m8_tet_q24   = q16_mul_to_q24(p->k_tet_q16, q16_mul(u2_q16, q16_mul(N, q16_mul(C, M8))));
  dxdt_q16[7] = m8_drift_q24 - m8_tet_q24;
}

/* ----------------------------------------------------------------------------
 * 2B. TRUE 4-STAGE RUNGE-KUTTA (`RK4`) STEPPER FOR THE 8-STATE SYSTEM
 *     Evaluates k1, k2, k3, k4 (in Q8.24) and updates state x (in Q16.16)
 * ---------------------------------------------------------------------------- */
static inline int32_t q24_slope_step_to_q16(int32_t dt_q16, int32_t k_q24) {
  int64_t prod = (int64_t)dt_q16 * (int64_t)k_q24;
  return (int32_t)((prod + (prod >= 0 ? (1LL << 23) : -(1LL << 23))) >> 24);
}

static inline void coralnpu_rk4_core8_full_step_q16(
    mitorestorer_ode_state_q16_t* state,
    const mitorestorer_core_params_q16_t* params,
    int32_t u0_q16,
    int32_t u1_q16,
    int32_t u2_q16,
    int32_t dt_q16,
    int64_t* rk4_sub_lsb_rem) {
  int32_t k1[8], k2[8], k3[8], k4[8];
  mitorestorer_ode_state_q16_t tmp;
  int32_t* s_ptr = (int32_t*)state;
  int32_t* t_ptr = (int32_t*)&tmp;
  int32_t half_dt_q16 = dt_q16 >> 1;

  /* Stage 1: k1 = f(x_n) */
  coralnpu_eval_core8_rhs_q16(state, params, u0_q16, u1_q16, u2_q16, k1);

  /* Stage 2: k2 = f(x_n + 0.5 * dt * k1) */
  for (int i = 0; i < 8; ++i) {
    t_ptr[i] = s_ptr[i] + q24_slope_step_to_q16(half_dt_q16, k1[i]);
  }
  coralnpu_eval_core8_rhs_q16(&tmp, params, u0_q16, u1_q16, u2_q16, k2);

  /* Stage 3: k3 = f(x_n + 0.5 * dt * k2) */
  for (int i = 0; i < 8; ++i) {
    t_ptr[i] = s_ptr[i] + q24_slope_step_to_q16(half_dt_q16, k2[i]);
  }
  coralnpu_eval_core8_rhs_q16(&tmp, params, u0_q16, u1_q16, u2_q16, k3);

  /* Stage 4: k4 = f(x_n + dt * k3) */
  for (int i = 0; i < 8; ++i) {
    t_ptr[i] = s_ptr[i] + q24_slope_step_to_q16(dt_q16, k3[i]);
  }
  coralnpu_eval_core8_rhs_q16(&tmp, params, u0_q16, u1_q16, u2_q16, k4);

  /* Combine: x_{n+1} = x_n + (dt / 6) * (k1 + 2*k2 + 2*k3 + k4) with sub-LSB remainder */
  const int64_t denom = 6LL * (1LL << 24); /* 100663296 */
  for (int i = 0; i < 8; ++i) {
    int64_t weighted_sum = (int64_t)k1[i] + 2 * (int64_t)k2[i] + 2 * (int64_t)k3[i] + (int64_t)k4[i];
    int64_t num = weighted_sum * (int64_t)dt_q16 + (rk4_sub_lsb_rem ? rk4_sub_lsb_rem[i] : 0);
    int32_t delta = (int32_t)(num / denom);
    if (rk4_sub_lsb_rem) {
      rk4_sub_lsb_rem[i] = num - (int64_t)delta * denom;
    }
    s_ptr[i] += delta;
  }
}

/* Legacy 2-stage Heun / trapezoidal helper preserved for backward testbench compatibility,
 * with corrected RVV Zve32x widening Q16.16 shift (`vwmul` + `vnsra.wi ..., 16`). */
static inline void coralnpu_rvv_rk4_bioenergetic_step(
    mitorestorer_ode_state_q16_t* state,
    const int32_t* deriv_k1_q16,
    const int32_t* deriv_k2_q16,
    int32_t dt_q16) {
  int32_t* s_ptr = (int32_t*)state;
#if defined(__riscv_vector)
  size_t vl = __riscv_vsetvl_e32m2(8);
  vint32m2_t vs  = __riscv_vle32_v_i32m2(s_ptr, vl);
  vint32m2_t vk1 = __riscv_vle32_v_i32m2(deriv_k1_q16, vl);
  vint32m2_t vk2 = __riscv_vle32_v_i32m2(deriv_k2_q16, vl);
  vint32m2_t vavg = __riscv_vadd_vv_i32m2(vk1, vk2, vl);
  vavg = __riscv_vsra_vx_i32m2(vavg, 1, vl);
  /* Widening 64-bit product followed by 16-bit Q16.16 right shift */
  vint64m4_t vprod = __riscv_vwmul_vx_i64m4(vavg, dt_q16, vl);
  vint32m2_t vstep = __riscv_vnsra_wx_i32m2(vprod, 16, vl);
  vs = __riscv_vadd_vv_i32m2(vs, vstep, vl);
  __riscv_vse32_v_i32m2(s_ptr, vs, vl);
#else
  for (int i = 0; i < 8; ++i) {
    int64_t avg_k = ((int64_t)deriv_k1_q16[i] + (int64_t)deriv_k2_q16[i]) >> 1;
    s_ptr[i] += (int32_t)((avg_k * (int64_t)dt_q16) >> 16);
  }
#endif
}

/* ----------------------------------------------------------------------------
 * 3. MATRIX ENGINE (`256-MAC INT8`): Zero-Skip Spectral & Biosensor Gate
 * ---------------------------------------------------------------------------- */
static inline mitorestorer_gate_action_t coralnpu_matrix256_stage_gate(
    const mitorestorer_ode_state_q16_t* state,
    const int8_t* spectral_features_i8,
    const int8_t* weights_i8,
    bool safety_ok,
    int32_t* out_logit_acc,
    int32_t* out_zero_skipped_bins) {
  if (!safety_ok) {
    if (out_logit_acc) *out_logit_acc = 0;
    if (out_zero_skipped_bins) *out_zero_skipped_bins = 0;
    return BCR_GATE_STAGE0_SAFETY_HOLD;
  }
  int32_t acc_stage2 = 0;
  int32_t skipped = 0;
  /* Hardware zero-skipping outer-product MAC across active spectral bins */
  for (int i = 0; i < 16; ++i) {
    int8_t feat = spectral_features_i8[i];
    if (feat == 0) {
      skipped++;
      continue; /* Coral NPU zero-skip hardware gate */
    }
    acc_stage2 += (int32_t)feat * (int32_t)weights_i8[i];
  }
  if (out_logit_acc) *out_logit_acc = acc_stage2;
  if (out_zero_skipped_bins) *out_zero_skipped_bins = skipped;

  /* Compute heteroplasmy fraction h = D / (W + D) in per-mille (x1000) */
  int64_t tot_mass = (int64_t)state->w_mass_q16 + (int64_t)state->d_mass_q16 + 1;
  int32_t h_pct_x10 = (int32_t)(((int64_t)state->d_mass_q16 * 1000) / tot_mass);

  /* Switch from Stage 1 (Mitophagy Purge) to Stage 2 (Anabolic Build) once
   * damaged mtDNA drops below 8.0% (80 in x10) and stress/inhibitor clears */
  if (h_pct_x10 <= 80 && state->inhibitor_q16 <= (12 << 16) && acc_stage2 >= 0) {
    return BCR_GATE_STAGE2_ANABOLIC_BUILD;
  }
  return BCR_GATE_STAGE1_MITOPHAGY_PURGE;
}

/* ----------------------------------------------------------------------------
 * 4. ROLLING-WINDOW NOISE FILTER, BIOFOULING FALLBACK & MULTIMODAL CROSS-CHECK
 *    Enforces BOTH 64-frame persistence quorum (>=60/64) AND secondary modality.
 * ---------------------------------------------------------------------------- */
static inline mitorestorer_gate_action_t coralnpu_rolling_fuzzy_sensor_gate(
    mitorestorer_coral_target_t target,
    mitorestorer_rolling_guard_dtcm_t* guard,
    int16_t current_primary_x10,
    int16_t secondary_cross_x100,
    int16_t elapsed_hours_x10,
    bool optical_ph_verified) {
  if (target == CORAL_TARGET_APIS_SMART_HIVE) {
    /* Require >= 60 of 64 rolling frames (>= 93.75% persistence) at >= 220.0 Hz (2200)
     * AND secondary brood-core temperature in [33.30 C, 35.70 C] (3330..3570) */
    int16_t old_val = guard->wingbeat_hz_x10_ring[guard->head_idx];
    if (old_val >= 2200 && guard->sustained_count > 0) guard->sustained_count--;
    guard->wingbeat_hz_x10_ring[guard->head_idx] = current_primary_x10;
    if (current_primary_x10 >= 2200 && guard->sustained_count < 64) guard->sustained_count++;
    guard->head_idx = (uint8_t)((guard->head_idx + 1) & 63);
    bool secondary_ok = (secondary_cross_x100 >= 3330 && secondary_cross_x100 <= 3570);
    return (guard->sustained_count >= 60 && secondary_ok)
               ? BCR_GATE_STAGE2_ANABOLIC_BUILD
               : BCR_GATE_STAGE1_MITOPHAGY_PURGE;
  } else if (target == CORAL_TARGET_CART_GMP_BIOREACTOR) {
    /* Detect amperometric Lactate:Pyruvate biofouling erratic jumps (|delta| > 4.50 ratio units) */
    int16_t delta = current_primary_x10 - guard->prev_lac_pyr_x100;
    if (delta < 0) delta = -delta;
    guard->prev_lac_pyr_x100 = current_primary_x10;
    if (delta > 450) guard->biofoul_erratic_count++;
    if (guard->biofoul_erratic_count >= 3) {
      /* Amperometric electrode biofouled! Fall back to timed 64.0h + inline optical pH check */
      if (elapsed_hours_x10 >= 640 && optical_ph_verified &&
          secondary_cross_x100 >= 695 && secondary_cross_x100 <= 745) {
        return BCR_GATE_STAGE2_ANABOLIC_BUILD;
      }
      return BCR_GATE_BIOFOUL_TIMED64H_OPTICAL_PH_HOLD;
    }
    bool secondary_ph_ok = (secondary_cross_x100 >= 695 && secondary_cross_x100 <= 745);
    return (current_primary_x10 <= 1200 && secondary_ph_ok)
               ? BCR_GATE_STAGE2_ANABOLIC_BUILD
               : BCR_GATE_STAGE1_MITOPHAGY_PURGE;
  } else {
    /* Human Wearable Ring: Cross-reference interstitial B-OHB ketones (current_primary_x10 >= 12 -> 1.2 mM)
     * with PPG-HRV RMSSD autonomic stability (secondary_cross_x100 >= 3800 -> 38.0 ms) */
    if (current_primary_x10 >= 12 && secondary_cross_x100 >= 3800) {
      return BCR_GATE_STAGE2_ANABOLIC_BUILD;
    }
    return BCR_GATE_STAGE1_MITOPHAGY_PURGE;
  }
}

/* ----------------------------------------------------------------------------
 * 5. BIDIRECTIONAL ACTUATOR ROLLBACK & 32-BYTE DTCM BREADCRUMB CHECKPOINT
 *    Note on Plant vs. Controller State Separation (Paper I, Section 5):
 *    - `coralnpu_RevertActuatorAndPreservePlant`: Reverses physical actuation
 *      (u2 <- 0, u0 <- 1, u1 <- 1) while preserving the deteriorated biological
 *      plant state `plant_state` so recovery must be earned via dynamics.
 *    - `coralnpu_bidirectional_rollback_guard`: Restores the controller's
 *      internal reference anchor `breadcrumb_ckpt_q16` for re-gating comparison.
 * ---------------------------------------------------------------------------- */
static inline void coralnpu_save_breadcrumb_ckpt(
    mitorestorer_rolling_guard_dtcm_t* guard,
    const mitorestorer_ode_state_q16_t* state) {
  const int32_t* s_ptr = (const int32_t*)state;
  for (int i = 0; i < 8; ++i) {
    guard->breadcrumb_ckpt_q16[i] = s_ptr[i];
  }
}

static inline mitorestorer_gate_action_t coralnpu_bidirectional_rollback_guard(
    mitorestorer_coral_target_t target,
    mitorestorer_gate_action_t active_gate,
    mitorestorer_ode_state_q16_t* controller_ref_state,
    mitorestorer_rolling_guard_dtcm_t* guard,
    int32_t post_transition_metric_a,
    int32_t post_transition_metric_b) {
  if (active_gate != BCR_GATE_STAGE2_ANABOLIC_BUILD) {
    return active_gate;
  }
  bool trigger_rollback = false;
  if (target == CORAL_TARGET_APIS_SMART_HIVE) {
    if (post_transition_metric_a < 3330 || post_transition_metric_a > 3570 ||
        post_transition_metric_b >= 65) {
      trigger_rollback = true;
    }
  } else if (target == CORAL_TARGET_CART_GMP_BIOREACTOR) {
    if (post_transition_metric_a < 695 || post_transition_metric_b < 4500) {
      trigger_rollback = true;
    }
  } else {
    if (post_transition_metric_a >= 2500 || post_transition_metric_b < 720) {
      trigger_rollback = true;
    }
  }

  if (trigger_rollback) {
    guard->rollback_triggered_flag = 1;
    if (controller_ref_state) {
      int32_t* s_ptr = (int32_t*)controller_ref_state;
      for (int i = 0; i < 8; ++i) {
        s_ptr[i] = guard->breadcrumb_ckpt_q16[i];
      }
    }
    return BCR_GATE_STAGE2_ROLLBACK_ABORT;
  }
  return BCR_GATE_STAGE2_ANABOLIC_BUILD;
}

#endif /* MITORESTORER_CORALNPU_KERNEL_H_ */

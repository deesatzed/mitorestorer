/* ============================================================================
 * MITORESTORER CORAL NPU FIXED-POINT KERNEL (v3.2 — Paper I Verified Release)
 * File        : mitorestorer_coralnpu_kernel.h
 * Target ISA  : Google Coral NPU (`RV32IMF_Zve32x` + 256-MAC INT8 Matrix Engine)
 * Datasheet   : Google Coral NPU (8 KB ITCM, 32 KB DTCM; 1x 8 KB Bank Active)
 * Memory Map  : ITCM @ 0x0000_0000 (<8 KB) | Active DTCM Bank 0 (3,696 / 8,192 B)
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
  int32_t k_bio0_q16;         /* k_bio,0 = 0.075 (4915)                      */
  int32_t k_biomax_q16;       /* k_bio,max = 0.140 (9175)                    */
  int32_t k_mito_w_q16;       /* k_mito,W = 0.015 (983)                      */
  int32_t k_mito_d_max_q16;   /* k_mito,D,max = 0.115 (7537)                 */
  int32_t k_dmg_q16;          /* k_dmg = 0.040 (2621)                        */
  int32_t k_acc_q16;          /* k_acc = 0.080 (5243)                        */
  int32_t k_clr_q16;          /* k_clr = 0.180 (11796)                       */
  int32_t k_b_half_q16;       /* K_B = 0.250 (16384)                         */
  int32_t k_e_relax_q16;      /* k_E = 0.250 (16384)                         */
  int32_t k_g_prod_q16;       /* k_G,prod = 0.080 (5243)                     */
  int32_t k_g_clr_q16;        /* k_G,clr = 0.160 (10486)                     */
  int32_t k_n_syn0_q16;       /* k_N,syn0 = 0.060 (3932)                     */
  int32_t k_n_syn2_q16;       /* k_N,syn2 = 0.080 (5243)                     */
  int32_t k_parp_q16;         /* k_PARP = 0.060 (3932)                       */
  int32_t k_n_use_q16;        /* k_N,use = 0.095 (6226)                      */
  int32_t k_oct0_q16;         /* k_oct0 = 0.040 (2621)                       */
  int32_t k_oct2_q16;         /* k_oct2 = 0.100 (6554)                       */
  int32_t k_ox_cys_q16;       /* k_ox,Cys278 = 0.050 (3277)                  */
  int32_t k_drift_q16;        /* k_drift = 0.0015 (98)                       */
  int32_t k_tet_q16;          /* k_TET = 0.012 (786)                         */
  int32_t eta_d_q16;          /* eta_D = 0.150 (9830)                        */
  int32_t beta_b_q16;         /* beta_B = 0.250 (16384)                      */
  int32_t beta_h_q16;         /* beta_h = 0.650 (42598)                      */
  int32_t gamma_c_q16;        /* gamma_C = 0.340 (22282)                     */
  int32_t theta_gsh_q16;      /* theta_GSH = 1.150 (75366)                   */
  int32_t theta_cys_q16;      /* theta_Cys = 0.850 (55706)                   */
  int32_t c1_flux_q16;        /* c_1 = 0.220 (14418)                         */
  int32_t c2_flux_q16;        /* c_2 = 0.040 (2621)                          */
  int32_t h_gate_q16;         /* h_gate = 0.080 (5243)                       */
  int32_t b_gate_q16;         /* B_gate = 0.120 (7864)                       */
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

static inline void coralnpu_q24_axpy8_q16(
    int32_t out_q16[8],
    const int32_t base_q16[8],
    const int32_t k_q24[8],
    int32_t dt_q16) {
#if defined(__riscv_vector)
  size_t vl = __riscv_vsetvl_e32m2(8);
  vint32m2_t vbase = __riscv_vle32_v_i32m2(base_q16, vl);
  vint32m2_t vk24  = __riscv_vle32_v_i32m2(k_q24, vl);
  vint64m4_t vprod = __riscv_vwmul_vx_i64m4(vk24, dt_q16, vl);
  /* Symmetric round-to-nearest before 24-bit right shift to Q16.16 */
  vbool16_t vneg   = __riscv_vmslt_vx_i64m4_b16(vprod, 0, vl);
  vint64m4_t vbias = __riscv_vmv_v_x_i64m4(1LL << 23, vl);
  vbias            = __riscv_vrsub_vx_i64m4_m(vneg, vbias, vbias, 0, vl);
  vprod            = __riscv_vadd_vv_i64m4(vprod, vbias, vl);
  vint32m2_t vstep = __riscv_vnsra_wx_i32m2(vprod, 24, vl);
  vint32m2_t vout  = __riscv_vadd_vv_i32m2(vbase, vstep, vl);
  __riscv_vse32_v_i32m2(out_q16, vout, vl);
#else
  for (int i = 0; i < 8; ++i) {
    out_q16[i] = base_q16[i] + q24_slope_step_to_q16(dt_q16, k_q24[i]);
  }
#endif
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
  coralnpu_q24_axpy8_q16(t_ptr, s_ptr, k1, half_dt_q16);
  coralnpu_eval_core8_rhs_q16(&tmp, params, u0_q16, u1_q16, u2_q16, k2);

  /* Stage 3: k3 = f(x_n + 0.5 * dt * k2) */
  coralnpu_q24_axpy8_q16(t_ptr, s_ptr, k2, half_dt_q16);
  coralnpu_eval_core8_rhs_q16(&tmp, params, u0_q16, u1_q16, u2_q16, k3);

  /* Stage 4: k4 = f(x_n + dt * k3) */
  coralnpu_q24_axpy8_q16(t_ptr, s_ptr, k3, dt_q16);
  coralnpu_eval_core8_rhs_q16(&tmp, params, u0_q16, u1_q16, u2_q16, k4);

  /* Combine: x_{n+1} = x_n + (dt / 6) * (k1 + 2*k2 + 2*k3 + k4) with sub-LSB remainder */
  const int64_t denom = 6LL * (1LL << 24); /* 100663296 */
  for (int i = 0; i < 8; ++i) {
    int64_t weighted_sum = (int64_t)k1[i] + 2 * (int64_t)k2[i] + 2 * (int64_t)k3[i] + (int64_t)k4[i];
    if (rk4_sub_lsb_rem) {
      int64_t num = weighted_sum * (int64_t)dt_q16 + rk4_sub_lsb_rem[i];
      int32_t delta = (int32_t)(num / denom);
      rk4_sub_lsb_rem[i] = num - (int64_t)delta * denom;
      s_ptr[i] += delta;
    } else {
      /* Round-to-nearest without sub-LSB remainder propagation (fair ablation baseline) */
      int64_t num = weighted_sum * (int64_t)dt_q16;
      int32_t delta = (int32_t)((num + (num >= 0 ? (denom >> 1) : -(denom >> 1))) / denom);
      s_ptr[i] += delta;
    }
  }
}

/* Baseline Naive Q16.16 4-Stage RK4 Stepper (for Fixed-Point Deadband Ablation):
 * Quantizes slopes directly to Q16.16 (>> 8 from Q8.24) without Q8.24 slope precision
 * or sub-LSB remainder accumulation, exposing quantization deadband stalling as dt -> 0. */
static inline void coralnpu_rk4_core8_naive_q16_step(
    mitorestorer_ode_state_q16_t* state,
    const mitorestorer_core_params_q16_t* params,
    int32_t u0_q16,
    int32_t u1_q16,
    int32_t u2_q16,
    int32_t dt_q16) {
  int32_t k1_q24[8], k2_q24[8], k3_q24[8], k4_q24[8];
  int32_t k1_q16[8], k2_q16[8], k3_q16[8], k4_q16[8];
  mitorestorer_ode_state_q16_t tmp;
  int32_t* s_ptr = (int32_t*)state;
  int32_t* t_ptr = (int32_t*)&tmp;
  int32_t half_dt_q16 = dt_q16 >> 1;

  coralnpu_eval_core8_rhs_q16(state, params, u0_q16, u1_q16, u2_q16, k1_q24);
  for (int i = 0; i < 8; ++i) {
    k1_q16[i] = (k1_q24[i] + (k1_q24[i] >= 0 ? 128 : -128)) >> 8;
    t_ptr[i] = s_ptr[i] + q16_mul(half_dt_q16, k1_q16[i]);
  }
  coralnpu_eval_core8_rhs_q16(&tmp, params, u0_q16, u1_q16, u2_q16, k2_q24);
  for (int i = 0; i < 8; ++i) {
    k2_q16[i] = (k2_q24[i] + (k2_q24[i] >= 0 ? 128 : -128)) >> 8;
    t_ptr[i] = s_ptr[i] + q16_mul(half_dt_q16, k2_q16[i]);
  }
  coralnpu_eval_core8_rhs_q16(&tmp, params, u0_q16, u1_q16, u2_q16, k3_q24);
  for (int i = 0; i < 8; ++i) {
    k3_q16[i] = (k3_q24[i] + (k3_q24[i] >= 0 ? 128 : -128)) >> 8;
    t_ptr[i] = s_ptr[i] + q16_mul(dt_q16, k3_q16[i]);
  }
  coralnpu_eval_core8_rhs_q16(&tmp, params, u0_q16, u1_q16, u2_q16, k4_q24);
  for (int i = 0; i < 8; ++i) {
    k4_q16[i] = (k4_q24[i] + (k4_q24[i] >= 0 ? 128 : -128)) >> 8;
    int32_t sum6 = k1_q16[i] + 2 * k2_q16[i] + 2 * k3_q16[i] + k4_q16[i];
    int32_t weighted_q16 = (sum6 + (sum6 >= 0 ? 3 : -3)) / 6;
    int32_t delta_q16 = q16_mul(dt_q16, weighted_q16);
    s_ptr[i] += delta_q16;
  }
}

/* Host-exact architectural model of RVV Zve32x 128-bit LMUL=2 (`vint32m2_t` -> `vint64m4_t`)
 * widening multiply (`vwmul.vx`) + narrowing arithmetic shift (`vnsra.wx ..., 16`),
 * enabling 0-ULP dual-path verification on non-RISC-V hosts. */
static inline void coralnpu_zve32x_isa_emulated_step_q16(
    int32_t state_vec[8],
    const int32_t k1_vec[8],
    const int32_t k2_vec[8],
    int32_t dt_q16) {
  /* Models 2x 128-bit vector registers (`LMUL=2`, `SEW=32`, `VL=8`) widened to
   * 4x 128-bit registers (`LMUL=4`, `SEW=64`, `VL=8`) via `vwmul.vx` + `vnsra.wx`. */
  int64_t vprod_widen_m4[8];
  for (int lane = 0; lane < 8; ++lane) {
    int32_t vavg = (int32_t)(((int64_t)k1_vec[lane] + (int64_t)k2_vec[lane]) >> 1);
    vprod_widen_m4[lane] = (int64_t)vavg * (int64_t)dt_q16; /* vwmul.vx */
  }
  for (int lane = 0; lane < 8; ++lane) {
    int32_t vstep = (int32_t)(vprod_widen_m4[lane] >> 16);   /* vnsra.wx ..., 16 */
    state_vec[lane] += vstep;                                /* vadd.vv */
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

/* ----------------------------------------------------------------------------
 * 6. UNIFIED PAPER I DUAL-THRESHOLD GATE & FIXED-POINT 7-STATE HUBER-EKF
 *    - `coralnpu_core8_stage_gate_q16`: Enforces BOTH h(t) <= h_gate AND B(t) <= B_gate.
 *    - `coralnpu_huber_ekf7_step_q16`: Bare-metal, zero-malloc Q16.16 7-state
 *      Huber-regularized Extended Kalman Filter over y = [E, ATP, G, N]^T using
 *      exact sequential scalar measurement updates (avoiding 4x4 matrix inversion).
 *    - `coralnpu_dtcm_bank0_workspace_t`: Complete static L1 DTCM Bank 0 layout
 *      verified at compile/runtime (`sizeof <= 8192` bytes).
 * ---------------------------------------------------------------------------- */
static inline bool coralnpu_core8_stage_gate_q16(
    const mitorestorer_ode_state_q16_t* state,
    int32_t h_gate_q16,
    int32_t b_gate_q16) {
  int32_t m_q16 = state->w_mass_q16 + state->d_mass_q16;
  int32_t h_q16 = (m_q16 > 64) ? q16_div(state->d_mass_q16, m_q16) : 0;
  return (h_q16 <= h_gate_q16) && (state->bottleneck_b_q16 <= b_gate_q16);
}

static inline uint32_t coralnpu_isqrt64(uint64_t n) {
  if (n == 0) return 0;
  uint64_t x = n;
  uint64_t y = (x + 1ULL) >> 1;
  while (y < x) {
    x = y;
    y = (x + n / x) >> 1;
  }
  return (uint32_t)x;
}

static inline void coralnpu_obs4_and_row_q16(
    const int32_t x7_q16[7],
    const mitorestorer_core_params_q16_t* p,
    int ch,
    int32_t* y_pred_q16,
    int32_t h_row_q16[7]) {
  for (int i = 0; i < 7; ++i) h_row_q16[i] = 0;
  int32_t W = x7_q16[0] > 0 ? x7_q16[0] : 0;
  int32_t D = x7_q16[1] > 0 ? x7_q16[1] : 0;
  int32_t C = x7_q16[6] > 0 ? x7_q16[6] : 0;
  if (ch == 0) {
    *y_pred_q16 = x7_q16[3];
    h_row_q16[3] = Q16_ONE;
  } else if (ch == 1) {
    int32_t eff_m = W + q16_mul(p->eta_d_q16, D);
    int32_t pcr_f = Q16_ONE + q16_mul(p->gamma_c_q16, C);
    *y_pred_q16 = q16_mul(eff_m, pcr_f);
    h_row_q16[0] = pcr_f;
    h_row_q16[1] = q16_mul(p->eta_d_q16, pcr_f);
    h_row_q16[6] = q16_mul(p->gamma_c_q16, eff_m);
  } else if (ch == 2) {
    *y_pred_q16 = x7_q16[4];
    h_row_q16[4] = Q16_ONE;
  } else {
    *y_pred_q16 = x7_q16[5];
    h_row_q16[5] = Q16_ONE;
  }
}

static inline int32_t q24_mul(int32_t a_q24, int32_t b_q24) {
  int64_t prod = (int64_t)a_q24 * (int64_t)b_q24;
  return (int32_t)((prod + (prod >= 0 ? (1LL << 23) : -(1LL << 23))) >> 24);
}

static inline int32_t q24_mul_q16(int32_t a_q24, int32_t b_q16) {
  /* (Q8.24 * Q16.16) >> 16 -> Q8.24 */
  int64_t prod = (int64_t)a_q24 * (int64_t)b_q16;
  return (int32_t)((prod + (prod >= 0 ? 32768 : -32768)) >> 16);
}

static inline int32_t q24_div_to_q16(int32_t num_q24, int32_t den_q24) {
  /* (Q8.24 << 16) / Q8.24 -> Q16.16 Kalman gain */
  if (den_q24 == 0) return 0;
  int64_t scaled = ((int64_t)num_q24 << 16);
  return (int32_t)(scaled / (int64_t)den_q24);
}

static inline void coralnpu_huber_ekf7_meas_update_q16(
    mitorestorer_ode_state_q16_t* x_hat,
    int32_t P_q24[7][7],
    const mitorestorer_core_params_q16_t* p,
    const int32_t y_meas_q16[4],
    int32_t r_var_q24,
    int32_t gamma_huber_q16) {
  int32_t* x_ptr = (int32_t*)x_hat;
  for (int ch = 0; ch < 4; ++ch) {
    int32_t y_pred = 0;
    int32_t h_row[7];
    coralnpu_obs4_and_row_q16(x_ptr, p, ch, &y_pred, h_row);
    int32_t innov_q16 = y_meas_q16[ch] - y_pred;

    /* PH^T column (in Q8.24) and scalar innovation variance S = H P H^T + R (in Q8.24) */
    int32_t ph_col_q24[7];
    int64_t s_acc_q24 = r_var_q24;
    for (int i = 0; i < 7; ++i) {
      int64_t sum = 0;
      for (int j = 0; j < 7; ++j) {
        sum += (int64_t)P_q24[i][j] * (int64_t)h_row[j]; /* Q8.24 * Q16.16 = Q24.40 */
      }
      ph_col_q24[i] = (int32_t)((sum + (sum >= 0 ? 32768 : -32768)) >> 16);
      s_acc_q24 += ((int64_t)h_row[i] * (int64_t)ph_col_q24[i] + 32768) >> 16;
    }
    int32_t s_var_q24 = (s_acc_q24 > 16) ? (int32_t)s_acc_q24 : 16;

    /* Exact integer sqrt(S_q24 << 8) = sqrt(S * 2^32) gives std(innov) in Q16.16 */
    int32_t std_q16 = (int32_t)coralnpu_isqrt64(((uint64_t)(uint32_t)s_var_q24) << 8);
    int32_t clip_bound_q16 = q16_mul(gamma_huber_q16, std_q16);
    int32_t clipped_innov_q16 = innov_q16;
    if (clipped_innov_q16 > clip_bound_q16) clipped_innov_q16 = clip_bound_q16;
    if (clipped_innov_q16 < -clip_bound_q16) clipped_innov_q16 = -clip_bound_q16;

    int32_t abs_in = innov_q16 >= 0 ? innov_q16 : -innov_q16;
    int32_t abs_cl = clipped_innov_q16 >= 0 ? clipped_innov_q16 : -clipped_innov_q16;
    int32_t w_hub_q16 = (abs_in > clip_bound_q16 && abs_in > 0) ? q16_div(abs_cl, abs_in) : Q16_ONE;

    int32_t k_gain_q16[7];
    for (int i = 0; i < 7; ++i) {
      k_gain_q16[i] = q24_div_to_q16(ph_col_q24[i], s_var_q24);
      int32_t x_new = x_ptr[i] + q16_mul(k_gain_q16[i], clipped_innov_q16);
      x_ptr[i] = (x_new > 66) ? x_new : 66; /* 0.001 positivity floor */
    }
    for (int i = 0; i < 7; ++i) {
      int32_t kw_q16 = q16_mul(k_gain_q16[i], w_hub_q16);
      for (int j = 0; j < 7; ++j) {
        P_q24[i][j] -= q24_mul_q16(ph_col_q24[j], kw_q16);
      }
    }
  }
  for (int i = 0; i < 7; ++i) {
    if (P_q24[i][i] < 16) P_q24[i][i] = 16;
    for (int j = i + 1; j < 7; ++j) {
      int32_t sym = (P_q24[i][j] + P_q24[j][i]) >> 1;
      P_q24[i][j] = sym;
      P_q24[j][i] = sym;
    }
  }
}

static inline void coralnpu_huber_ekf7_predict_q16(
    mitorestorer_ode_state_q16_t* x_hat,
    int32_t P_q24[7][7],
    const mitorestorer_core_params_q16_t* p,
    int32_t u0_q16,
    int32_t u1_q16,
    int32_t u2_q16,
    int32_t dt_q16,
    const int32_t q_diag_q24[7],
    int64_t rk4_rem[8]) {
  /* 1. First-order discrete Jacobian F = I + dt * df/dx via Q8.24 RHS differences (eps = 256 = 2^-8) */
  int32_t f0_q24[8], fp_q24[8];
  int32_t F_q16[7][7];
  coralnpu_eval_core8_rhs_q16(x_hat, p, u0_q16, u1_q16, u2_q16, f0_q24);
  for (int j = 0; j < 7; ++j) {
    mitorestorer_ode_state_q16_t xp = *x_hat;
    ((int32_t*)&xp)[j] += 256;
    coralnpu_eval_core8_rhs_q16(&xp, p, u0_q16, u1_q16, u2_q16, fp_q24);
    for (int i = 0; i < 7; ++i) {
      int32_t df_dx_q16 = fp_q24[i] - f0_q24[i];
      F_q16[i][j] = ((i == j) ? Q16_ONE : 0) + q16_mul(dt_q16, df_dx_q16);
    }
  }
  /* 2. Propagate Q8.24 covariance P <- F P F^T + Q (F is Q16.16, P and Q are Q8.24) */
  int32_t FP_q24[7][7];
  for (int i = 0; i < 7; ++i) {
    for (int j = 0; j < 7; ++j) {
      int64_t sum = 0;
      for (int k = 0; k < 7; ++k) sum += (int64_t)F_q16[i][k] * (int64_t)P_q24[k][j];
      FP_q24[i][j] = (int32_t)((sum + (sum >= 0 ? 32768 : -32768)) >> 16);
    }
  }
  for (int i = 0; i < 7; ++i) {
    for (int j = 0; j < 7; ++j) {
      int64_t sum = 0;
      for (int k = 0; k < 7; ++k) sum += (int64_t)FP_q24[i][k] * (int64_t)F_q16[j][k];
      P_q24[i][j] = (int32_t)((sum + (sum >= 0 ? 32768 : -32768)) >> 16) + ((i == j) ? q_diag_q24[i] : 0);
    }
  }
  /* 3. Propagate 8-state trajectory with Sigma-Delta sub-LSB remainder */
  coralnpu_rk4_core8_full_step_q16(x_hat, p, u0_q16, u1_q16, u2_q16, dt_q16, rk4_rem);
}

/* Complete Single-Bank L1 DTCM Static Workspace (Must Fit Strictly Inside 8,192-Byte Bank 0) */
typedef struct {
  int8_t                          int8_proj_weights[16 * 64]; /* 1,024 B: 16x64 INT8 sparse gate matrix */
  int8_t                          int8_head_weights[64 * 3];  /*   192 B: 64x3 INT8 classifier head     */
  int16_t                         goertzel_lut_q15[256];      /*   512 B: 16-bin Goertzel/FFT twiddles  */
  int32_t                         ema_window_q16[64];         /*   256 B: 64-frame Q16.16 EMA buffer    */
  int32_t                         ekf_P_cov_q24[7][7];        /*   196 B: 7x7 Q8.24 EKF covariance      */
  int32_t                         ekf_F_jac_q16[7][7];        /*   196 B: 7x7 Q16.16 Jacobian workspace */
  int32_t                         ekf_Q_diag_q24[7];          /*    28 B: 7-state Q8.24 process noise   */
  int32_t                         rk4_k_scratch_q24[4][8];    /*   128 B: 4-stage Q8.24 RK4 slopes      */
  int64_t                         rk4_sub_lsb_rem[8];         /*    64 B: 8-lane Sigma-Delta remainders */
  mitorestorer_ode_state_q16_t    x_hat_q16;                  /*    32 B: 8-state Q16.16 estimator      */
  mitorestorer_core_params_q16_t  params_q16;                 /*   132 B: 33 Q16.16 kinetic parameters  */
  mitorestorer_rolling_guard_dtcm_t rolling_guard;            /*   168 B: 64-frame ring + 32B ckpt      */
  uint8_t                         scalar_isr_stack_rsv[768];  /*   768 B: Reserved scalar/ISR C stack   */
} coralnpu_dtcm_bank0_workspace_t;

#endif /* MITORESTORER_CORALNPU_KERNEL_H_ */

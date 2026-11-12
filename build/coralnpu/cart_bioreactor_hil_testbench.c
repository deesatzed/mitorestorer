/* ============================================================================
 * Host-Emulated C99 Fixed-Point CAR-T Perfusion Cassette Verification Testbench
 * File: build/coralnpu/cart_bioreactor_hil_testbench.c
 * Executes:
 *   1. Exact 1-to-1 8-state CAR-T bioreactor ODE parity: IEEE-754 FP64 vs.
 *      Scalar FP32 vs. Q16.16 + Q8.24 slopes + 64-bit Sigma-Delta RK4 vs.
 *      Naive Q16.16 RK4 across N in {500, 1000, 2000, 4000, 8000} steps.
 *   2. 100% Integer FPU-Free Sequential 5-Channel (R_LP, OUR_vol, OD, pH, ORR)
 *      Huber-EKF + NIS > 13.5 Biofouling Handover vs. IEEE-754 FP64 EKF under
 *      donor kinetic mismatch and progressive amperometric electrode fouling.
 *   3. Post-switch cassette perfusion occlusion acidosis rescue via autonomous
 *      neutral bypass flush (u_flush = 1) with continuous plant state x(t).
 *   4. Zero-heap static SRAM workspace verification (< 8,192 B, no unused NPU
 *      neural weight padding).
 * ============================================================================ */

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "mitorestorer_coralnpu_kernel.h"

/* Q12.20 Fixed-Point Helpers for EKF Covariance Matrix P (Range [-2048, +2048), LSB = 2^-20) */
#define Q20_FROM_FLOAT(x) ((int32_t)((x) >= 0.0 ? ((x) * 1048576.0 + 0.5) : ((x) * 1048576.0 - 0.5)))

static inline uint32_t isqrt64_u32(uint64_t n) {
  uint64_t x = n;
  uint64_t y = (x + 1ULL) >> 1;
  while (y < x) {
    x = y;
    y = (x + n / x) >> 1;
  }
  return (x > 0xFFFFFFFFULL) ? 0xFFFFFFFFU : (uint32_t)x;
}

/* 8-State CAR-T Perfusion Cassette Parameter Struct in Q16.16 */
typedef struct {
  int32_t K_q16;
  int32_t k_mito_w_q16;
  int32_t k_bio_0_q16;
  int32_t k_bio_max_q16;
  int32_t k_mito_d_max_q16;
  int32_t k_rep_adv_q16;
  int32_t xi_d_q16;
  int32_t k_dmg_q16;
  int32_t theta_gsh_q16;
  int32_t R_basal_q16;
  int32_t K_LP_q16;
  int32_t eta_d_q16;
  int32_t gamma_pdh_q16;
  int32_t gamma_src_q16;
  int32_t alpha_gal_q16;
  int32_t kappa_glyc_q16;
  int32_t c_akt_q16;
  int32_t K_R_ATP_q16;
  int32_t beta_lp_q16;
  int32_t beta_h_q16;
  int32_t beta_ph_q16;
  int32_t c1_q16;
  int32_t c2_q16;
  int32_t mu_glyc_q16;
  int32_t mu_ox_q16;
  int32_t K_ATP_q16;
  int32_t K_W_div_q16;
  int32_t N_cap_q16;
  int32_t k_warburg_q16;
  int32_t alpha_h_q16;
  int32_t k_pdh_q16;
  int32_t beta_azd_q16;
  int32_t k_perf_lp_q16;
  int32_t k_flush_lp_q16;
  int32_t R_media_q16;
  int32_t k_imp_q16;
  int32_t f_0_q16;
  int32_t alpha_akt_q16;
  int32_t phi_sirt1_q16;
  int32_t k_exp_q16;
  int32_t gamma_acid_q16;
  int32_t s_basal_q16;
  int32_t k_cristae_q16;
  int32_t beta_src_q16;
  int32_t k_src_dec_q16;
  int32_t k_exh_q16;
  int32_t theta_exh_q16;
  int32_t k_warb_exh_q16;
  int32_t k_acid_exh_q16;
  int32_t r_epig_0_q16;
  int32_t k_rev_epig_q16;
  int32_t k_perf_0_q16;
  int32_t k_perf_flush_q16;
  int32_t pH_media_q16;
  int32_t k_acid_q16;
  int32_t R_buf_q16;
  int32_t N_ref_q16;
  int32_t K_SRC_mem_q16;
} cart_cassette_params_q16_t;

/* Deterministic Zero-Heap Static SRAM CAR-T Bioreactor Workspace (< 8,192 B) */
typedef struct {
  int16_t                         optical_lut_q15[256];       /* 512 B: Dual-wavelength pH/optode LUT  */
  int32_t                         innov_window_q16[64];       /* 256 B: 64-frame innovation ring       */
  int32_t                         ekf_P_cov_q20[8][8];        /* 256 B: 8x8 Q12.20 EKF covariance      */
  int32_t                         ekf_F_jac_q16[8][8];        /* 256 B: 8x8 Q16.16 Jacobian scratch    */
  int32_t                         ekf_Q_diag_q20[8];          /*  32 B: 8-state Q12.20 process noise   */
  int32_t                         ekf_R_var_q20[5];           /*  20 B: 5-ch Q12.20 base sensor var    */
  int32_t                         n_at_24h_q16;               /*   4 B: 24h OD reference fold          */
  uint8_t                         nis_streak;                 /*   1 B: Consecutive NIS > 13.5 count   */
  uint8_t                         ch0_disabled;               /*   1 B: Biofouled R_LP channel latch   */
  uint8_t                         stage;                      /*   1 B: Active cassette stage (1 or 2) */
  uint8_t                         quorum;                     /*   1 B: 4-frame persistence counter    */
  uint32_t                        reserved_align;             /*   4 B: 8-byte alignment               */
  int32_t                         rk4_k_scratch_q24[4][8];    /* 128 B: 4-stage Q8.24 RK4 slopes       */
  int64_t                         rk4_sub_lsb_rem[8];         /*  64 B: 8-lane Sigma-Delta remainder   */
  int32_t                         x_hat_q16[8];               /*  32 B: 8-state Q16.16 EKF estimator   */
  cart_cassette_params_q16_t      params_q16;                 /* 232 B: 58 Q16.16 cassette params      */
  mitorestorer_rolling_guard_dtcm_t rolling_guard;            /* 168 B: Persistence guard & anchor     */
  uint8_t                         scalar_isr_stack_rsv[736];  /* 736 B: Reserved deterministic C stack */
} coralnpu_cart_dtcm_bank0_workspace_t;

static void init_cart_params_q16(cart_cassette_params_q16_t* p) {
  p->K_q16            = Q16_FROM_FLOAT(1.40);
  p->k_mito_w_q16     = Q16_FROM_FLOAT(0.008);
  p->k_bio_0_q16      = Q16_FROM_FLOAT(0.028);
  p->k_bio_max_q16    = Q16_FROM_FLOAT(0.086);
  p->k_mito_d_max_q16 = Q16_FROM_FLOAT(0.068);
  p->k_rep_adv_q16    = Q16_FROM_FLOAT(1.22);
  p->xi_d_q16         = Q16_FROM_FLOAT(1.00);
  p->k_dmg_q16        = Q16_FROM_FLOAT(0.0012);
  p->theta_gsh_q16    = Q16_FROM_FLOAT(1.38);
  p->R_basal_q16      = Q16_FROM_FLOAT(8.0);
  p->K_LP_q16         = Q16_FROM_FLOAT(22.0);
  p->eta_d_q16        = Q16_FROM_FLOAT(0.22);
  p->gamma_pdh_q16    = Q16_FROM_FLOAT(0.72);
  p->gamma_src_q16    = Q16_FROM_FLOAT(0.22);
  p->alpha_gal_q16    = Q16_FROM_FLOAT(0.74);
  p->kappa_glyc_q16   = Q16_FROM_FLOAT(1.40);
  p->c_akt_q16        = Q16_FROM_FLOAT(0.08);
  p->K_R_ATP_q16      = Q16_FROM_FLOAT(8.0);
  p->beta_lp_q16      = Q16_FROM_FLOAT(0.40);
  p->beta_h_q16       = Q16_FROM_FLOAT(0.88);
  p->beta_ph_q16      = Q16_FROM_FLOAT(1.90);
  p->c1_q16           = Q16_FROM_FLOAT(0.24);
  p->c2_q16           = Q16_FROM_FLOAT(0.06);
  p->mu_glyc_q16      = Q16_FROM_FLOAT(0.062);
  p->mu_ox_q16        = Q16_FROM_FLOAT(0.0185);
  p->K_ATP_q16        = Q16_FROM_FLOAT(0.52);
  p->K_W_div_q16      = Q16_FROM_FLOAT(0.56);
  p->N_cap_q16        = Q16_FROM_FLOAT(85.0);
  p->k_warburg_q16    = Q16_FROM_FLOAT(3.00);
  p->alpha_h_q16      = Q16_FROM_FLOAT(0.46);
  p->k_pdh_q16        = Q16_FROM_FLOAT(0.088);
  p->beta_azd_q16     = Q16_FROM_FLOAT(2.80);
  p->k_perf_lp_q16    = Q16_FROM_FLOAT(0.120);
  p->k_flush_lp_q16   = Q16_FROM_FLOAT(0.85);
  p->R_media_q16      = Q16_FROM_FLOAT(5.0);
  p->k_imp_q16        = Q16_FROM_FLOAT(0.056);
  p->f_0_q16          = Q16_FROM_FLOAT(0.30);
  p->alpha_akt_q16    = Q16_FROM_FLOAT(0.30);
  p->phi_sirt1_q16    = Q16_FROM_FLOAT(1.18);
  p->k_exp_q16        = Q16_FROM_FLOAT(0.032);
  p->gamma_acid_q16   = Q16_FROM_FLOAT(7.50);
  p->s_basal_q16      = Q16_FROM_FLOAT(0.044);
  p->k_cristae_q16    = Q16_FROM_FLOAT(0.074);
  p->beta_src_q16     = Q16_FROM_FLOAT(1.32);
  p->k_src_dec_q16    = Q16_FROM_FLOAT(0.048);
  p->k_exh_q16        = Q16_FROM_FLOAT(0.0038);
  p->theta_exh_q16    = Q16_FROM_FLOAT(1.36);
  p->k_warb_exh_q16   = Q16_FROM_FLOAT(0.0032);
  p->k_acid_exh_q16   = Q16_FROM_FLOAT(0.220);
  p->r_epig_0_q16     = Q16_FROM_FLOAT(0.0032);
  p->k_rev_epig_q16   = Q16_FROM_FLOAT(0.0110);
  p->k_perf_0_q16     = Q16_FROM_FLOAT(0.48);
  p->k_perf_flush_q16 = Q16_FROM_FLOAT(1.80);
  p->pH_media_q16     = Q16_FROM_FLOAT(7.36);
  p->k_acid_q16       = Q16_FROM_FLOAT(0.0070);
  p->R_buf_q16        = Q16_FROM_FLOAT(12.0);
  p->N_ref_q16        = Q16_FROM_FLOAT(48.0);
  p->K_SRC_mem_q16    = Q16_FROM_FLOAT(0.68);
}

/* IEEE-754 FP64 Reference 8-State CAR-T Cassette RHS (With Optional Donor Multipliers) */
static void eval_cart8_rhs_fp64_ex(
    const double x[8],
    double u1,
    double u2,
    double u_fl,
    bool fault_occl,
    double mu_glyc_mult,
    double k_mito_d_mult,
    double k_warb_mult,
    double dxdt[8]) {
  const double K = 1.40, k_mito_w = 0.008, k_bio_0 = 0.028, k_bio_max = 0.086;
  const double k_mito_d_max = 0.068 * k_mito_d_mult, k_rep_adv = 1.22, xi_d = 1.00;
  const double k_dmg = 0.0012, theta_gsh = 1.38, R_basal = 8.0, K_LP = 22.0;
  const double eta_d = 0.22, gamma_pdh = 0.72, gamma_src = 0.22, alpha_gal = 0.74;
  const double kappa_glyc = 1.40, c_akt = 0.08, K_R_ATP = 8.0;
  const double beta_lp = 0.40, beta_h = 0.88, beta_ph = 1.90, c1 = 0.24, c2 = 0.06;
  const double mu_glyc = 0.062 * mu_glyc_mult, mu_ox = 0.0185, K_ATP = 0.52, K_W_div = 0.56, N_cap = 85.0;
  const double k_warburg = 3.00 * k_warb_mult, alpha_h = 0.46, k_pdh = 0.088, beta_azd = 2.80;
  const double k_perf_lp = 0.120, k_flush_lp = 0.85, R_media = 5.0;
  const double k_imp = 0.056, f_0 = 0.30, alpha_akt = 0.30, phi_sirt1 = 1.18;
  const double k_exp = 0.032, gamma_acid = 7.50;
  const double s_basal = 0.044, k_cristae = 0.074, beta_src = 1.32, k_src_dec = 0.048;
  const double k_exh = 0.0038, theta_exh = 1.36, k_warb_exh = 0.0032, k_acid_exh = 0.220;
  const double r_epig_0 = 0.0032, k_rev_epig = 0.0110;
  const double k_perf_0 = 0.48, k_perf_flush = 1.80, pH_media = 7.36;
  const double k_acid = 0.0070, R_buf = 12.0, N_ref = 48.0;

  double W = x[0] > 0.0 ? x[0] : 0.0;
  double D = x[1] > 0.0 ? x[1] : 0.0;
  double N = x[2] > 0.0 ? x[2] : 0.0;
  double R = x[3] > 0.0 ? x[3] : 0.0;
  double F = x[4] > 0.0 ? (x[4] < 1.0 ? x[4] : 1.0) : 0.0;
  double S = x[5] > 0.0 ? x[5] : 0.0;
  double E = x[6] > 0.0 ? (x[6] < 1.0 ? x[6] : 1.0) : 0.0;
  double pH = x[7];

  double M = (W + D > 1e-9) ? (W + D) : 1e-9;
  double h = D / M;
  double w_frac = W / M;

  double g_glc = fault_occl ? 1.0 : (1.0 - alpha_gal * u2);
  double u_azd_eff = fault_occl ? 0.0 : u2;
  double u_sirt = u2;

  double r_exc = (R > R_basal) ? (R - R_basal) : 0.0;
  double v_unmask = K_LP / (K_LP + r_exc);
  double k_md = k_mito_w + (k_mito_d_max - k_mito_w) * u1 * v_unmask;

  double s_exc = (S > 1.0) ? (S - 1.0) : 0.0;
  double atp_ox = (W + eta_d * D) * (1.0 + gamma_pdh * u_azd_eff) * (1.0 + gamma_src * s_exc);
  double atp_glyc = kappa_glyc * g_glc * (1.0 - c_akt * u1 * (1.0 - u2)) * (R / (R + K_R_ATP));
  double atp_tot = atp_ox + atp_glyc;

  double ph_def = (7.02 > pH) ? (7.02 - pH) : 0.0;
  double ratio_ph = ph_def / 0.16;
  double phi_ph = 1.0 / (1.0 + ratio_ph * ratio_ph);

  double eff_m = W + eta_d * D;
  if (eff_m < 0.26) eff_m = 0.26;
  double e_dot = (1.0 + beta_lp * (r_exc / K_LP) + beta_h * h + beta_ph * ph_def) / eff_m;
  double f_eff = 1.0 - c1 * u1 + c2 * u2;
  double ep = e_dot * f_eff * f_eff;
  double dmg = (ep > theta_gsh) ? (k_dmg * (ep - theta_gsh)) : 0.0;

  double crowd = (M < K) ? (1.0 - M / K) : 0.0;
  dxdt[0] = (k_bio_0 + (k_bio_max - k_bio_0) * u_sirt * F) * W * crowd - (k_mito_w + dmg) * W;
  dxdt[1] = k_rep_adv * (k_bio_0 + xi_d * (k_bio_max - k_bio_0) * u2) * D * crowd - k_md * D + dmg * W;

  double w_func = W * (1.0 - h);
  double ox_comp = (w_func * w_func) / (w_func * w_func + K_W_div * K_W_div);
  double pdh_gate = 0.15 + 0.85 * u_azd_eff;
  double mu_anab = mu_glyc * g_glc * (1.0 - 0.08 * u1 * (1.0 - u2)) + mu_ox * u2 * pdh_gate * ox_comp;
  double n_cap_term = (N < N_cap) ? (1.0 - N / N_cap) : 0.0;
  dxdt[2] = mu_anab * phi_ph * (atp_tot / (atp_tot + K_ATP)) * N * n_cap_term;

  double perf_lp_eff = fault_occl ? (0.08 * k_perf_lp) : k_perf_lp;
  dxdt[3] = k_warburg * g_glc * (1.0 - 0.35 * u1 * (1.0 - u2)) * (1.0 + alpha_h * h)
          - k_pdh * (1.0 + beta_azd * u_azd_eff) * w_frac * R
          - perf_lp_eff * (R - R_media)
          - k_flush_lp * u_fl * (R - R_media);

  double acid_f = (6.90 > pH) ? (6.90 - pH) : 0.0;
  dxdt[4] = k_imp * (f_0 + alpha_akt * u1 + phi_sirt1 * u_sirt * pdh_gate) * phi_ph * (1.0 - F)
          - k_exp * (1.0 - 0.52 * u1) * (1.0 - 0.85 * u_sirt * pdh_gate) * (1.0 + gamma_acid * acid_f) * F;

  dxdt[5] = (s_basal * (W + eta_d * D) + k_cristae * w_func * (1.0 + beta_src * u_sirt * F * pdh_gate))
          - k_src_dec * S;

  double ep_exh = (ep > theta_exh) ? (ep - theta_exh) : 0.0;
  double acid_exh = (6.86 > pH) ? (6.86 - pH) : 0.0;
  dxdt[6] = (k_exh * ep_exh + k_warb_exh * g_glc * (1.0 - 0.22 * u1) * (1.0 + 0.70 * (N / N_ref)) + k_acid_exh * acid_exh) * (1.0 - E)
          - (r_epig_0 + k_rev_epig * u_sirt * F * pdh_gate) * E;

  double k_perf_eff = (fault_occl ? (0.08 * k_perf_0) : k_perf_0) + k_perf_flush * u_fl;
  double r_buf_exc = (R > R_buf) ? (R - R_buf) : 0.0;
  double acid_mult = fault_occl ? 1.65 : g_glc;
  dxdt[7] = k_perf_eff * (pH_media - pH) - k_acid * acid_mult * r_buf_exc * (N / N_ref);
}

static void eval_cart8_rhs_fp64(const double x[8], double u1, double u2, double u_fl, bool fault_occl, double dxdt[8]) {
  eval_cart8_rhs_fp64_ex(x, u1, u2, u_fl, fault_occl, 1.0, 1.0, 1.0, dxdt);
}

static void rk4_cart8_step_fp64_ex(
    double x[8], double u1, double u2, double u_fl, bool fault, double dt,
    double mu_g_m, double k_md_m, double k_wb_m) {
  double k1[8], k2[8], k3[8], k4[8], tmp[8];
  eval_cart8_rhs_fp64_ex(x, u1, u2, u_fl, fault, mu_g_m, k_md_m, k_wb_m, k1);
  for (int i = 0; i < 8; ++i) tmp[i] = x[i] + 0.5 * dt * k1[i];
  eval_cart8_rhs_fp64_ex(tmp, u1, u2, u_fl, fault, mu_g_m, k_md_m, k_wb_m, k2);
  for (int i = 0; i < 8; ++i) tmp[i] = x[i] + 0.5 * dt * k2[i];
  eval_cart8_rhs_fp64_ex(tmp, u1, u2, u_fl, fault, mu_g_m, k_md_m, k_wb_m, k3);
  for (int i = 0; i < 8; ++i) tmp[i] = x[i] + dt * k3[i];
  eval_cart8_rhs_fp64_ex(tmp, u1, u2, u_fl, fault, mu_g_m, k_md_m, k_wb_m, k4);
  for (int i = 0; i < 8; ++i) {
    x[i] += (dt / 6.0) * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]);
  }
}

static void rk4_cart8_step_fp64(double x[8], double u1, double u2, double u_fl, bool fault, double dt) {
  rk4_cart8_step_fp64_ex(x, u1, u2, u_fl, fault, dt, 1.0, 1.0, 1.0);
}

/* Scalar IEEE-754 FP32 Reference Stepper */
static void eval_cart8_rhs_fp32(const float x[8], float u1, float u2, float u_fl, bool fault_occl, float dxdt[8]) {
  double xd[8], dxdtd[8];
  for (int i = 0; i < 8; ++i) xd[i] = (double)x[i];
  eval_cart8_rhs_fp64(xd, (double)u1, (double)u2, (double)u_fl, fault_occl, dxdtd);
  for (int i = 0; i < 8; ++i) dxdt[i] = (float)dxdtd[i];
}

static void rk4_cart8_step_fp32(float x[8], float u1, float u2, float u_fl, bool fault, float dt) {
  float k1[8], k2[8], k3[8], k4[8], tmp[8];
  eval_cart8_rhs_fp32(x, u1, u2, u_fl, fault, k1);
  for (int i = 0; i < 8; ++i) tmp[i] = x[i] + 0.5f * dt * k1[i];
  eval_cart8_rhs_fp32(tmp, u1, u2, u_fl, fault, k2);
  for (int i = 0; i < 8; ++i) tmp[i] = x[i] + 0.5f * dt * k2[i];
  eval_cart8_rhs_fp32(tmp, u1, u2, u_fl, fault, k3);
  for (int i = 0; i < 8; ++i) tmp[i] = x[i] + dt * k3[i];
  eval_cart8_rhs_fp32(tmp, u1, u2, u_fl, fault, k4);
  for (int i = 0; i < 8; ++i) {
    x[i] += (dt / 6.0f) * (k1[i] + 2.0f * k2[i] + 2.0f * k3[i] + k4[i]);
  }
}

/* Bare-Metal Q16.16 State -> Q8.24 Slope Evaluator (100% FPU-Free Integer Math) */
static void eval_cart8_rhs_q16_to_q24(
    const int32_t x_q16[8],
    const cart_cassette_params_q16_t* p,
    int32_t u1_q16,
    int32_t u2_q16,
    int32_t u_fl_q16,
    bool fault_occl,
    int32_t dxdt_q24[8]) {
  int32_t W = q16_max0(x_q16[0]);
  int32_t D = q16_max0(x_q16[1]);
  int32_t N = q16_max0(x_q16[2]);
  int32_t R = q16_max0(x_q16[3]);
  int32_t F = x_q16[4] < 0 ? 0 : (x_q16[4] > Q16_ONE ? Q16_ONE : x_q16[4]);
  int32_t S = q16_max0(x_q16[5]);
  int32_t E = x_q16[6] < 0 ? 0 : (x_q16[6] > Q16_ONE ? Q16_ONE : x_q16[6]);
  int32_t pH = x_q16[7];

  int32_t M = W + D;
  if (M < 64) M = 64;
  int32_t h = q16_div(D, M);
  int32_t w_frac = q16_div(W, M);

  int32_t g_glc = fault_occl ? Q16_ONE : (Q16_ONE - q16_mul(p->alpha_gal_q16, u2_q16));
  int32_t u_azd_eff = fault_occl ? 0 : u2_q16;
  int32_t u_sirt = u2_q16;

  int32_t r_exc = q16_max0(R - p->R_basal_q16);
  int32_t v_unmask = q16_div(p->K_LP_q16, p->K_LP_q16 + r_exc);
  int32_t k_md = p->k_mito_w_q16 + q16_mul(p->k_mito_d_max_q16 - p->k_mito_w_q16, q16_mul(u1_q16, v_unmask));

  int32_t eff_mass_raw = W + q16_mul(p->eta_d_q16, D);
  int32_t s_exc = q16_max0(S - Q16_ONE);
  int32_t atp_ox = q16_mul(
      eff_mass_raw,
      q16_mul(Q16_ONE + q16_mul(p->gamma_pdh_q16, u_azd_eff),
              Q16_ONE + q16_mul(p->gamma_src_q16, s_exc)));
  int32_t u1_only = q16_mul(u1_q16, Q16_ONE - u2_q16);
  int32_t atp_glyc = q16_mul(
      q16_mul(p->kappa_glyc_q16, g_glc),
      q16_mul(Q16_ONE - q16_mul(p->c_akt_q16, u1_only),
              q16_div(R, R + p->K_R_ATP_q16)));
  int32_t atp_tot = atp_ox + atp_glyc;

  int32_t ph_def = q16_max0(Q16_FROM_FLOAT(7.02) - pH);
  int32_t ratio_ph = q16_div(ph_def, Q16_FROM_FLOAT(0.16));
  int32_t phi_ph = q16_div(Q16_ONE, Q16_ONE + q16_mul(ratio_ph, ratio_ph));

  int32_t eff_m = eff_mass_raw < Q16_FROM_FLOAT(0.26) ? Q16_FROM_FLOAT(0.26) : eff_mass_raw;
  int32_t num_edot = Q16_ONE
                   + q16_mul(p->beta_lp_q16, q16_div(r_exc, p->K_LP_q16))
                   + q16_mul(p->beta_h_q16, h)
                   + q16_mul(p->beta_ph_q16, ph_def);
  int32_t e_dot = q16_div(num_edot, eff_m);
  int32_t f_eff = Q16_ONE - q16_mul(p->c1_q16, u1_q16) + q16_mul(p->c2_q16, u2_q16);
  int32_t ep = q16_mul(e_dot, q16_mul(f_eff, f_eff));
  int32_t ep_exc = q16_max0(ep - p->theta_gsh_q16);
  int32_t dmg = q16_mul(p->k_dmg_q16, ep_exc);

  int32_t crowd = q16_max0(Q16_ONE - q16_div(M, p->K_q16));
  int32_t k_bio_w = p->k_bio_0_q16 + q16_mul(p->k_bio_max_q16 - p->k_bio_0_q16, q16_mul(u_sirt, F));
  dxdt_q24[0] = q16_mul_to_q24(k_bio_w, q16_mul(W, crowd)) - q16_mul_to_q24(p->k_mito_w_q16 + dmg, W);

  int32_t k_bio_d = q16_mul(p->k_rep_adv_q16, p->k_bio_0_q16 + q16_mul(p->xi_d_q16, q16_mul(p->k_bio_max_q16 - p->k_bio_0_q16, u2_q16)));
  dxdt_q24[1] = q16_mul_to_q24(k_bio_d, q16_mul(D, crowd)) - q16_mul_to_q24(k_md, D) + q16_mul_to_q24(dmg, W);

  int32_t w_func = q16_mul(W, Q16_ONE - h);
  int32_t wf2 = q16_mul(w_func, w_func);
  int32_t kw2 = q16_mul(p->K_W_div_q16, p->K_W_div_q16);
  int32_t ox_comp = q16_div(wf2, wf2 + kw2);
  int32_t pdh_gate = Q16_FROM_FLOAT(0.15) + q16_mul(Q16_FROM_FLOAT(0.85), u_azd_eff);
  int32_t mu_anab = q16_mul(p->mu_glyc_q16, q16_mul(g_glc, Q16_ONE - q16_mul(Q16_FROM_FLOAT(0.08), u1_only)))
                  + q16_mul(p->mu_ox_q16, q16_mul(u2_q16, q16_mul(pdh_gate, ox_comp)));
  int32_t atp_gate = q16_div(atp_tot, atp_tot + p->K_ATP_q16);
  int32_t n_crowd = q16_max0(Q16_ONE - q16_div(N, p->N_cap_q16));
  int32_t mu_eff = q16_mul(mu_anab, q16_mul(phi_ph, atp_gate));
  dxdt_q24[2] = q16_mul_to_q24(mu_eff, q16_mul(N, n_crowd));

  int32_t warb_gen = q16_mul(
      p->k_warburg_q16,
      q16_mul(g_glc, q16_mul(Q16_ONE - q16_mul(Q16_FROM_FLOAT(0.35), u1_only),
                             Q16_ONE + q16_mul(p->alpha_h_q16, h))));
  int32_t pdh_clr = q16_mul(p->k_pdh_q16, q16_mul(Q16_ONE + q16_mul(p->beta_azd_q16, u_azd_eff), w_frac));
  int32_t perf_lp = fault_occl ? q16_mul(Q16_FROM_FLOAT(0.08), p->k_perf_lp_q16) : p->k_perf_lp_q16;
  dxdt_q24[3] = (warb_gen << 8)
              - q16_mul_to_q24(pdh_clr, R)
              - q16_mul_to_q24(perf_lp + q16_mul(p->k_flush_lp_q16, u_fl_q16), R - p->R_media_q16);

  int32_t sirt_pdh = q16_mul(u_sirt, pdh_gate);
  int32_t imp_drv = p->f_0_q16 + q16_mul(p->alpha_akt_q16, u1_q16) + q16_mul(p->phi_sirt1_q16, sirt_pdh);
  int32_t acid_f = q16_max0(Q16_FROM_FLOAT(6.90) - pH);
  int32_t exp_drv = q16_mul(
      p->k_exp_q16,
      q16_mul(Q16_ONE - q16_mul(Q16_FROM_FLOAT(0.52), u1_q16),
              q16_mul(Q16_ONE - q16_mul(Q16_FROM_FLOAT(0.85), sirt_pdh),
                      Q16_ONE + q16_mul(p->gamma_acid_q16, acid_f))));
  dxdt_q24[4] = q16_mul_to_q24(q16_mul(p->k_imp_q16, q16_mul(imp_drv, phi_ph)), Q16_ONE - F)
              - q16_mul_to_q24(exp_drv, F);

  int32_t s_gen = q16_mul(p->s_basal_q16, eff_mass_raw)
                + q16_mul(p->k_cristae_q16, q16_mul(w_func, Q16_ONE + q16_mul(p->beta_src_q16, q16_mul(F, sirt_pdh))));
  dxdt_q24[5] = (s_gen << 8) - q16_mul_to_q24(p->k_src_dec_q16, S);

  int32_t ep_exh = q16_max0(ep - p->theta_exh_q16);
  int32_t acid_exh = q16_max0(Q16_FROM_FLOAT(6.86) - pH);
  int32_t n_ratio = q16_div(N, p->N_ref_q16);
  int32_t warb_exh = q16_mul(
      p->k_warb_exh_q16,
      q16_mul(g_glc, q16_mul(Q16_ONE - q16_mul(Q16_FROM_FLOAT(0.22), u1_q16),
                             Q16_ONE + q16_mul(Q16_FROM_FLOAT(0.70), n_ratio))));
  int32_t exh_in = q16_mul(p->k_exh_q16, ep_exh) + warb_exh + q16_mul(p->k_acid_exh_q16, acid_exh);
  int32_t exh_rev = p->r_epig_0_q16 + q16_mul(p->k_rev_epig_q16, q16_mul(F, sirt_pdh));
  dxdt_q24[6] = q16_mul_to_q24(exh_in, Q16_ONE - E) - q16_mul_to_q24(exh_rev, E);

  int32_t k_perf_eff = (fault_occl ? q16_mul(Q16_FROM_FLOAT(0.08), p->k_perf_0_q16) : p->k_perf_0_q16)
                     + q16_mul(p->k_perf_flush_q16, u_fl_q16);
  int32_t r_buf_exc = q16_max0(R - p->R_buf_q16);
  int32_t acid_mult = fault_occl ? Q16_FROM_FLOAT(1.65) : g_glc;
  int32_t acid_drop = q16_mul(p->k_acid_q16, q16_mul(acid_mult, q16_mul(r_buf_exc, n_ratio)));
  dxdt_q24[7] = q16_mul_to_q24(k_perf_eff, p->pH_media_q16 - pH) - (acid_drop << 8);
}

/* 4-Stage RK4 with Q8.24 Slopes and 64-bit Sigma-Delta Sub-LSB Remainder Accumulator */
static void rk4_cart8_step_q16_sd(
    int32_t x_q16[8],
    const cart_cassette_params_q16_t* p,
    int32_t u1_q16,
    int32_t u2_q16,
    int32_t u_fl_q16,
    bool fault,
    int32_t dt_q16,
    int64_t rem[8]) {
  int32_t k1[8], k2[8], k3[8], k4[8], tmp[8];
  eval_cart8_rhs_q16_to_q24(x_q16, p, u1_q16, u2_q16, u_fl_q16, fault, k1);
  for (int i = 0; i < 8; ++i) {
    int64_t dx = ((int64_t)dt_q16 * (int64_t)k1[i]) >> 25;
    tmp[i] = x_q16[i] + (int32_t)dx;
  }
  eval_cart8_rhs_q16_to_q24(tmp, p, u1_q16, u2_q16, u_fl_q16, fault, k2);
  for (int i = 0; i < 8; ++i) {
    int64_t dx = ((int64_t)dt_q16 * (int64_t)k2[i]) >> 25;
    tmp[i] = x_q16[i] + (int32_t)dx;
  }
  eval_cart8_rhs_q16_to_q24(tmp, p, u1_q16, u2_q16, u_fl_q16, fault, k3);
  for (int i = 0; i < 8; ++i) {
    int64_t dx = ((int64_t)dt_q16 * (int64_t)k3[i]) >> 24;
    tmp[i] = x_q16[i] + (int32_t)dx;
  }
  eval_cart8_rhs_q16_to_q24(tmp, p, u1_q16, u2_q16, u_fl_q16, fault, k4);
  const int64_t denom = 6LL << 24;
  for (int i = 0; i < 8; ++i) {
    int64_t wsum_q24 = (int64_t)k1[i] + 2LL * (int64_t)k2[i] + 2LL * (int64_t)k3[i] + (int64_t)k4[i];
    int64_t num = (int64_t)dt_q16 * wsum_q24 + rem[i];
    int64_t dq = num / denom;
    rem[i] = num - dq * denom;
    x_q16[i] += (int32_t)dq;
  }
}

/* Naive Q16.16 4-Stage RK4 (Truncated Q16.16 Slopes Without Sigma-Delta Diffusion) */
static void rk4_cart8_step_q16_naive(
    int32_t x_q16[8],
    const cart_cassette_params_q16_t* p,
    int32_t u1_q16,
    int32_t u2_q16,
    int32_t u_fl_q16,
    bool fault,
    int32_t dt_q16) {
  int32_t k1_24[8], k2_24[8], k3_24[8], k4_24[8], k1[8], k2[8], k3[8], k4[8], tmp[8];
  eval_cart8_rhs_q16_to_q24(x_q16, p, u1_q16, u2_q16, u_fl_q16, fault, k1_24);
  for (int i = 0; i < 8; ++i) {
    k1[i] = k1_24[i] >> 8;
    tmp[i] = x_q16[i] + (q16_mul(dt_q16, k1[i]) >> 1);
  }
  eval_cart8_rhs_q16_to_q24(tmp, p, u1_q16, u2_q16, u_fl_q16, fault, k2_24);
  for (int i = 0; i < 8; ++i) {
    k2[i] = k2_24[i] >> 8;
    tmp[i] = x_q16[i] + (q16_mul(dt_q16, k2[i]) >> 1);
  }
  eval_cart8_rhs_q16_to_q24(tmp, p, u1_q16, u2_q16, u_fl_q16, fault, k3_24);
  for (int i = 0; i < 8; ++i) {
    k3[i] = k3_24[i] >> 8;
    tmp[i] = x_q16[i] + q16_mul(dt_q16, k3[i]);
  }
  eval_cart8_rhs_q16_to_q24(tmp, p, u1_q16, u2_q16, u_fl_q16, fault, k4_24);
  for (int i = 0; i < 8; ++i) {
    k4[i] = k4_24[i] >> 8;
    int32_t wavg = (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]) / 6;
    x_q16[i] += q16_mul(dt_q16, wavg);
  }
}

/* ============================================================================
 * 5-Channel Cassette Sensor Observation & Sequential Huber-EKF + NIS Handover
 * Channels:
 *   y[0]: Amperometric [Lac]/[Pyr] ratio R_LP (subject to protein biofouling)
 *   y[1]: Dissolved-O2 optode volumetric respiration OUR_vol = N_cells * ATP_ox
 *   y[2]: NIR optical backscatter cell-density fold OD = alpha_OD * N_cells
 *   y[3]: Dual-wavelength phenol-red pH
 *   y[4]: Quartz-window autofluorescence Optical Redox Ratio ORR (cell-density invariant)
 * ============================================================================ */
static void cassette_obs_fp64(const double x[8], double u2, double c_orr_glyc, double alpha_od, double y[5]) {
  double W = x[0] > 0.0 ? x[0] : 0.0;
  double D = x[1] > 0.0 ? x[1] : 0.0;
  double N = x[2] > 0.0 ? x[2] : 0.0;
  double R = x[3] > 0.0 ? x[3] : 0.0;
  double S = x[5] > 0.0 ? x[5] : 0.0;
  double pH = x[7];
  double s_exc = (S > 1.0) ? (S - 1.0) : 0.0;
  double atp_ox = (W + 0.22 * D) * (1.0 + 0.72 * u2) * (1.0 + 0.22 * s_exc);
  double our_vol = N * atp_ox;
  double g_glc = 1.0 - 0.74 * u2;
  double den = W + D + c_orr_glyc * g_glc;
  if (den < 1e-6) den = 1e-6;
  y[0] = R;
  y[1] = our_vol;
  y[2] = alpha_od * N;
  y[3] = pH;
  y[4] = W / den;
}

static void cassette_obs_q16(
    const int32_t x_q16[8],
    const cart_cassette_params_q16_t* p,
    int32_t u2_q16,
    int32_t y_q16[5]) {
  int32_t W = q16_max0(x_q16[0]);
  int32_t D = q16_max0(x_q16[1]);
  int32_t N = q16_max0(x_q16[2]);
  int32_t R = q16_max0(x_q16[3]);
  int32_t S = q16_max0(x_q16[5]);
  int32_t pH = x_q16[7];
  int32_t s_exc = q16_max0(S - Q16_ONE);
  int32_t eff_m = W + q16_mul(p->eta_d_q16, D);
  int32_t g_p = Q16_ONE + q16_mul(p->gamma_pdh_q16, u2_q16);
  int32_t g_s = Q16_ONE + q16_mul(p->gamma_src_q16, s_exc);
  int32_t atp_ox = q16_mul(eff_m, q16_mul(g_p, g_s));
  int32_t our_vol = q16_mul(N, atp_ox);
  int32_t g_glc = Q16_ONE - q16_mul(p->alpha_gal_q16, u2_q16);
  int32_t c_g = q16_mul(Q16_FROM_FLOAT(0.35), g_glc);
  int32_t den = W + D + c_g;
  if (den < 64) den = 64;
  y_q16[0] = R;
  y_q16[1] = our_vol;
  y_q16[2] = N;
  y_q16[3] = pH;
  y_q16[4] = q16_div(W, den);
}

/* 100% Integer FPU-Free Sequential 5-Channel Huber-EKF + NIS > 13.5 Biofouling Handover */
static void cart_ekf_step_q16(
    coralnpu_cart_dtcm_bank0_workspace_t* ws,
    const int32_t y_meas_q16[5],
    int32_t u1_q16,
    int32_t u2_q16,
    int32_t dt_q16,
    bool* newly_disabled_out) {
  *newly_disabled_out = false;
  int32_t x_prev[8];
  for (int i = 0; i < 8; ++i) x_prev[i] = ws->x_hat_q16[i];

  /* 1. Compute central-difference discrete Jacobian F = I + dt * df/dx in Q16.16 (delta_x = 32 LSB = 2^-11) */
  for (int j = 0; j < 8; ++j) {
    int32_t xp[8], xm[8], fp_q24[8], fm_q24[8];
    for (int i = 0; i < 8; ++i) { xp[i] = x_prev[i]; xm[i] = x_prev[i]; }
    xp[j] += 32;
    xm[j] -= 32;
    eval_cart8_rhs_q16_to_q24(xp, &ws->params_q16, u1_q16, u2_q16, 0, false, fp_q24);
    eval_cart8_rhs_q16_to_q24(xm, &ws->params_q16, u1_q16, u2_q16, 0, false, fm_q24);
    for (int i = 0; i < 8; ++i) {
      /* (fp - fm) in Q8.24 over 2*delta_x = 2^-10 -> shift << 2 yields exact Q16.16 df_i/dx_j */
      int32_t df_dx_q16 = (fp_q24[i] - fm_q24[i]) << 2;
      ws->ekf_F_jac_q16[i][j] = q16_mul(df_dx_q16, dt_q16) + ((i == j) ? Q16_ONE : 0);
    }
  }

  /* 2. Predict State via Q16.16 + Q8.24 + Sigma-Delta RK4 */
  rk4_cart8_step_q16_sd(ws->x_hat_q16, &ws->params_q16, u1_q16, u2_q16, 0, false, dt_q16, ws->rk4_sub_lsb_rem);

  /* 3. Predict Covariance P_pred = F * P * F^T + Q in Q12.20 */
  int32_t FP_q20[8][8];
  for (int i = 0; i < 8; ++i) {
    for (int j = 0; j < 8; ++j) {
      int64_t acc = 0;
      for (int k = 0; k < 8; ++k) {
        acc += (int64_t)ws->ekf_F_jac_q16[i][k] * (int64_t)ws->ekf_P_cov_q20[k][j];
      }
      FP_q20[i][j] = (int32_t)((acc + (acc >= 0 ? 32768LL : -32768LL)) >> 16);
    }
  }
  for (int i = 0; i < 8; ++i) {
    for (int j = 0; j < 8; ++j) {
      int64_t acc = 0;
      for (int k = 0; k < 8; ++k) {
        acc += (int64_t)FP_q20[i][k] * (int64_t)ws->ekf_F_jac_q16[j][k];
      }
      int32_t p_ij = (int32_t)((acc + (acc >= 0 ? 32768LL : -32768LL)) >> 16);
      if (i == j) p_ij += ws->ekf_Q_diag_q20[i];
      ws->ekf_P_cov_q20[i][j] = p_ij;
    }
  }

  /* 4. Sequential Scalar Huber-EKF Updates Across 5 Channels */
  const cart_cassette_params_q16_t* p = &ws->params_q16;
  for (int m_ch = 0; m_ch < 5; ++m_ch) {
    int32_t y_hat_q16[5];
    cassette_obs_q16(ws->x_hat_q16, p, u2_q16, y_hat_q16);
    int32_t innov_q16 = y_meas_q16[m_ch] - y_hat_q16[m_ch];

    int32_t H_q16[8] = {0, 0, 0, 0, 0, 0, 0, 0};
    int32_t r_var_q20 = ws->ekf_R_var_q20[m_ch];
    if (m_ch == 0) {
      H_q16[3] = Q16_ONE;
    } else if (m_ch == 1) {
      int32_t W = q16_max0(ws->x_hat_q16[0]);
      int32_t D = q16_max0(ws->x_hat_q16[1]);
      int32_t N = q16_max0(ws->x_hat_q16[2]);
      int32_t S = q16_max0(ws->x_hat_q16[5]);
      int32_t s_ex = q16_max0(S - Q16_ONE);
      int32_t g_p = Q16_ONE + q16_mul(p->gamma_pdh_q16, u2_q16);
      int32_t g_s = Q16_ONE + q16_mul(p->gamma_src_q16, s_ex);
      int32_t gp_gs = q16_mul(g_p, g_s);
      int32_t eff_m = W + q16_mul(p->eta_d_q16, D);
      H_q16[0] = q16_mul(N, gp_gs);
      H_q16[1] = q16_mul(N, q16_mul(p->eta_d_q16, gp_gs));
      H_q16[2] = q16_mul(eff_m, gp_gs);
      if (S > Q16_ONE) {
        H_q16[5] = q16_mul(N, q16_mul(eff_m, q16_mul(g_p, p->gamma_src_q16)));
      }
      /* Volumetric optode noise variance: 0.05^2 + (0.045 * N_hat)^2 in Q12.20 */
      int32_t sig_prop_q16 = q16_mul(Q16_FROM_FLOAT(0.045), N);
      int64_t prop_var_q20 = ((int64_t)sig_prop_q16 * (int64_t)sig_prop_q16) >> 12;
      r_var_q20 += (int32_t)prop_var_q20;
    } else if (m_ch == 2) {
      H_q16[2] = Q16_ONE;
      /* OD backscatter noise variance: 0.35^2 + (0.04 * N_hat)^2 in Q12.20 */
      int32_t N = q16_max0(ws->x_hat_q16[2]);
      int32_t sig_od_q16 = q16_mul(Q16_FROM_FLOAT(0.040), N);
      int64_t od_var_q20 = ((int64_t)sig_od_q16 * (int64_t)sig_od_q16) >> 12;
      r_var_q20 += (int32_t)od_var_q20;
    } else if (m_ch == 3) {
      H_q16[7] = Q16_ONE;
    } else {
      int32_t W = q16_max0(ws->x_hat_q16[0]);
      int32_t D = q16_max0(ws->x_hat_q16[1]);
      int32_t c_g = q16_mul(Q16_FROM_FLOAT(0.35), Q16_ONE - q16_mul(p->alpha_gal_q16, u2_q16));
      int32_t den = W + D + c_g;
      if (den < 64) den = 64;
      int32_t den2 = q16_mul(den, den);
      if (den2 < 16) den2 = 16;
      H_q16[0] = q16_div(D + c_g, den2);
      H_q16[1] = -q16_div(W, den2);
    }

    int32_t PHt_q20[8];
    int64_t h_pht_acc = 0;
    for (int i = 0; i < 8; ++i) {
      int64_t acc = 0;
      for (int j = 0; j < 8; ++j) {
        acc += (int64_t)ws->ekf_P_cov_q20[i][j] * (int64_t)H_q16[j];
      }
      PHt_q20[i] = (int32_t)((acc + (acc >= 0 ? 32768LL : -32768LL)) >> 16);
      h_pht_acc += (int64_t)H_q16[i] * (int64_t)PHt_q20[i];
    }
    int32_t S_inn_q20 = (int32_t)((h_pht_acc + (h_pht_acc >= 0 ? 32768LL : -32768LL)) >> 16) + r_var_q20;
    if (S_inn_q20 < 16) S_inn_q20 = 16;

    int64_t innov_sq_q32 = (int64_t)innov_q16 * (int64_t)innov_q16;
    int64_t s_inn_q32 = (int64_t)S_inn_q20 << 12;

    if (m_ch == 0) {
      /* Exact integer NIS > 13.5 (= 27/2) biofouling residual check */
      bool bad_nis = (2LL * innov_sq_q32 > 27LL * s_inn_q32);
      ws->nis_streak = bad_nis ? (uint8_t)(ws->nis_streak + 1U) : 0U;
      if (!ws->ch0_disabled && ws->nis_streak >= 3U) {
        ws->ch0_disabled = 1U;
        *newly_disabled_out = true;
      }
      if (ws->ch0_disabled) {
        continue; /* Channel 0 permanently disabled; optical quartet handles estimation */
      }
    }

    int32_t huber_w_q16 = Q16_ONE;
    if (4LL * innov_sq_q32 > 25LL * s_inn_q32) { /* |innov| > 2.5 * sqrt(S_inn) */
      uint32_t std_q16 = isqrt64_u32((uint64_t)s_inn_q32);
      int32_t abs_inn = innov_q16 >= 0 ? innov_q16 : -innov_q16;
      huber_w_q16 = q16_div(q16_mul(Q16_FROM_FLOAT(2.5), (int32_t)std_q16), abs_inn);
    }

    int32_t K_q16[8];
    for (int i = 0; i < 8; ++i) {
      int32_t k_raw_q16 = (int32_t)(((int64_t)PHt_q20[i] << 16) / (int64_t)S_inn_q20);
      K_q16[i] = q16_mul(k_raw_q16, huber_w_q16);
      ws->x_hat_q16[i] += q16_mul(K_q16[i], innov_q16);
    }
    for (int i = 0; i < 7; ++i) {
      if (ws->x_hat_q16[i] < 0) ws->x_hat_q16[i] = 0;
    }
    if (ws->x_hat_q16[4] > Q16_ONE) ws->x_hat_q16[4] = Q16_ONE;
    if (ws->x_hat_q16[6] > Q16_ONE) ws->x_hat_q16[6] = Q16_ONE;

    int32_t P_old_q20[8][8];
    for (int i = 0; i < 8; ++i) {
      for (int j = 0; j < 8; ++j) P_old_q20[i][j] = ws->ekf_P_cov_q20[i][j];
    }
    for (int i = 0; i < 8; ++i) {
      for (int j = 0; j < 8; ++j) {
        int64_t kh_p_acc = 0;
        for (int k = 0; k < 8; ++k) {
          int32_t kh_ik_q16 = q16_mul(K_q16[i], H_q16[k]);
          kh_p_acc += (int64_t)kh_ik_q16 * (int64_t)P_old_q20[k][j];
        }
        ws->ekf_P_cov_q20[i][j] = P_old_q20[i][j]
            - (int32_t)((kh_p_acc + (kh_p_acc >= 0 ? 32768LL : -32768LL)) >> 16);
      }
    }
  }
}

/* Reference IEEE-754 FP64 Sequential 5-Channel Huber-EKF + NIS > 13.5 Handover */
static void cart_ekf_step_fp64(
    double x_hat[8],
    double P_cov[8][8],
    const double Q_diag[8],
    const double R_base_var[5],
    const double y_meas[5],
    double u1,
    double u2,
    double dt,
    int* nis_streak,
    bool* ch0_disabled,
    bool* newly_disabled_out) {
  *newly_disabled_out = false;
  double F_jac[8][8];
  const double eps = 488.28125e-6; /* 32 / 65536 = 2^-11 */
  for (int j = 0; j < 8; ++j) {
    double xp[8], xm[8], fp[8], fm[8];
    for (int i = 0; i < 8; ++i) { xp[i] = x_hat[i]; xm[i] = x_hat[i]; }
    xp[j] += eps; xm[j] -= eps;
    eval_cart8_rhs_fp64(xp, u1, u2, 0.0, false, fp);
    eval_cart8_rhs_fp64(xm, u1, u2, 0.0, false, fm);
    for (int i = 0; i < 8; ++i) {
      F_jac[i][j] = ((fp[i] - fm[i]) / (2.0 * eps)) * dt + ((i == j) ? 1.0 : 0.0);
    }
  }
  rk4_cart8_step_fp64(x_hat, u1, u2, 0.0, false, dt);

  double FP[8][8], P_pred[8][8];
  for (int i = 0; i < 8; ++i) {
    for (int j = 0; j < 8; ++j) {
      double s = 0.0;
      for (int k = 0; k < 8; ++k) s += F_jac[i][k] * P_cov[k][j];
      FP[i][j] = s;
    }
  }
  for (int i = 0; i < 8; ++i) {
    for (int j = 0; j < 8; ++j) {
      double s = 0.0;
      for (int k = 0; k < 8; ++k) s += FP[i][k] * F_jac[j][k];
      P_pred[i][j] = s + ((i == j) ? Q_diag[i] : 0.0);
    }
  }
  for (int i = 0; i < 8; ++i) {
    for (int j = 0; j < 8; ++j) P_cov[i][j] = P_pred[i][j];
  }

  for (int m_ch = 0; m_ch < 5; ++m_ch) {
    double y_hat[5];
    cassette_obs_fp64(x_hat, u2, 0.35, 1.0, y_hat);
    double innov = y_meas[m_ch] - y_hat[m_ch];
    double H[8] = {0};
    double r_var = R_base_var[m_ch];
    if (m_ch == 0) {
      H[3] = 1.0;
    } else if (m_ch == 1) {
      double W = x_hat[0] > 0.0 ? x_hat[0] : 0.0;
      double D = x_hat[1] > 0.0 ? x_hat[1] : 0.0;
      double N = x_hat[2] > 0.0 ? x_hat[2] : 0.0;
      double S = x_hat[5] > 0.0 ? x_hat[5] : 0.0;
      double s_ex = S > 1.0 ? (S - 1.0) : 0.0;
      double g_p = 1.0 + 0.72 * u2;
      double g_s = 1.0 + 0.22 * s_ex;
      double eff_m = W + 0.22 * D;
      H[0] = N * g_p * g_s;
      H[1] = N * 0.22 * g_p * g_s;
      H[2] = eff_m * g_p * g_s;
      if (S > 1.0) H[5] = N * eff_m * g_p * 0.22;
      r_var += (0.045 * N) * (0.045 * N);
    } else if (m_ch == 2) {
      H[2] = 1.0;
      double N = x_hat[2] > 0.0 ? x_hat[2] : 0.0;
      r_var += (0.040 * N) * (0.040 * N);
    } else if (m_ch == 3) {
      H[7] = 1.0;
    } else {
      double W = x_hat[0] > 0.0 ? x_hat[0] : 0.0;
      double D = x_hat[1] > 0.0 ? x_hat[1] : 0.0;
      double c_g = 0.35 * (1.0 - 0.74 * u2);
      double den = W + D + c_g;
      if (den < 1e-6) den = 1e-6;
      H[0] = (D + c_g) / (den * den);
      H[1] = -W / (den * den);
    }

    double PHt[8], S_inn = r_var;
    for (int i = 0; i < 8; ++i) {
      double s = 0.0;
      for (int j = 0; j < 8; ++j) s += P_cov[i][j] * H[j];
      PHt[i] = s;
      S_inn += H[i] * s;
    }
    if (S_inn < 1e-9) S_inn = 1e-9;

    if (m_ch == 0) {
      double nis = (innov * innov) / S_inn;
      *nis_streak = (nis > 13.5) ? (*nis_streak + 1) : 0;
      if (!(*ch0_disabled) && *nis_streak >= 3) {
        *ch0_disabled = true;
        *newly_disabled_out = true;
      }
      if (*ch0_disabled) continue;
    }

    double std_inn = sqrt(S_inn);
    double hw = (fabs(innov) <= 2.5 * std_inn) ? 1.0 : (2.5 * std_inn / fabs(innov));
    double K[8];
    for (int i = 0; i < 8; ++i) {
      K[i] = (PHt[i] / S_inn) * hw;
      x_hat[i] += K[i] * innov;
    }
    for (int i = 0; i < 7; ++i) if (x_hat[i] < 0.0) x_hat[i] = 0.0;
    if (x_hat[4] > 1.0) x_hat[4] = 1.0;
    if (x_hat[6] > 1.0) x_hat[6] = 1.0;

    double P_old[8][8];
    for (int i = 0; i < 8; ++i) for (int j = 0; j < 8; ++j) P_old[i][j] = P_cov[i][j];
    for (int i = 0; i < 8; ++i) {
      for (int j = 0; j < 8; ++j) {
        double s = 0.0;
        for (int k = 0; k < 8; ++k) s += (K[i] * H[k]) * P_old[k][j];
        P_cov[i][j] = P_old[i][j] - s;
      }
    }
  }
}

/* Deterministic SplitMix64 + Box-Muller PRNG for Reproducible C99 Sensor Noise */
static uint64_t splitmix64_next(uint64_t* state) {
  uint64_t z = (*state += 0x9E3779B97F4A7C15ULL);
  z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
  z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
  return z ^ (z >> 31);
}

static double prng_uniform01(uint64_t* state) {
  return (double)(splitmix64_next(state) >> 11) * (1.0 / 9007199254740992.0);
}

static double prng_normal01(uint64_t* state) {
  double u1 = prng_uniform01(state);
  double u2 = prng_uniform01(state);
  if (u1 < 1e-12) u1 = 1e-12;
  return sqrt(-2.0 * log(u1)) * cos(6.28318530717958647692 * u2);
}

static double compute_tscm_pct(const double x[8]) {
  return 100.0 * x[4] * (1.0 - x[6]) * (x[5] / (x[5] + 0.68));
}

int main(int argc, char** argv) {
  bool json_out = (argc > 1 && strcmp(argv[1], "--json") == 0);

  static coralnpu_dtcm_bank0_workspace_t core_ws;
  static coralnpu_cart_dtcm_bank0_workspace_t cart_ws;
  memset(&core_ws, 0, sizeof(core_ws));
  memset(&cart_ws, 0, sizeof(cart_ws));
  init_cart_params_q16(&cart_ws.params_q16);

  size_t core_bytes = sizeof(core_ws);
  size_t cart_bytes = sizeof(cart_ws);
  assert(core_bytes == 3696);
  assert(cart_bytes < 8192);

  const int step_counts[5] = {500, 1000, 2000, 4000, 8000};
  double err_sd_pct[5], err_naive_pct[5], err_fp32_pct[5];
  double sw_fp64_1000 = -1.0, sw_q16_1000 = -1.0;
  double n_fp64_1000 = 0.0, n_q16_1000 = 0.0;
  double tscm_fp64_1000 = 0.0, tscm_q16_1000 = 0.0;
  double h_fp64_1000 = 0.0, h_q16_1000 = 0.0;

  for (int s = 0; s < 5; ++s) {
    int n_steps = step_counts[s];
    double dt = 168.0 / (double)n_steps;
    int32_t dt_q16 = Q16_FROM_FLOAT(dt);

    double x64[8] = {0.54, 0.42, 1.0, 22.0, 0.28, 1.30, 0.08, 7.32};
    float  x32[8] = {0.54f, 0.42f, 1.0f, 22.0f, 0.28f, 1.30f, 0.08f, 7.32f};
    int32_t x_sd[8], x_nv[8];
    int64_t rem[8] = {0};
    for (int i = 0; i < 8; ++i) {
      x_sd[i] = Q16_FROM_FLOAT(x64[i]);
      x_nv[i] = x_sd[i];
    }

    int st64 = 1, st32 = 1, st_sd = 1, st_nv = 1;
    double sw64 = -1.0, sw_sd = -1.0;

    for (int k = 0; k < n_steps; ++k) {
      double t = k * dt;
      double h64 = x64[1] / (x64[0] + x64[1]);
      if (st64 == 1 && h64 <= 0.045 && x64[2] >= 14.0 && x64[7] >= 7.10) {
        st64 = 2; sw64 = t;
      }
      double u1_64 = (st64 == 1) ? 1.0 : 0.0;
      double u2_64 = (st64 == 1) ? 0.0 : 1.0;
      rk4_cart8_step_fp64(x64, u1_64, u2_64, 0.0, false, dt);

      float h32 = x32[1] / (x32[0] + x32[1]);
      if (st32 == 1 && h32 <= 0.045f && x32[2] >= 14.0f && x32[7] >= 7.10f) {
        st32 = 2;
      }
      rk4_cart8_step_fp32(x32, st32 == 1 ? 1.0f : 0.0f, st32 == 1 ? 0.0f : 1.0f, 0.0f, false, (float)dt);

      int32_t h_sd = q16_div(x_sd[1], x_sd[0] + x_sd[1]);
      if (st_sd == 1 && h_sd <= Q16_FROM_FLOAT(0.045) && x_sd[2] >= Q16_FROM_FLOAT(14.0) && x_sd[7] >= Q16_FROM_FLOAT(7.10)) {
        st_sd = 2; sw_sd = t;
      }
      rk4_cart8_step_q16_sd(x_sd, &cart_ws.params_q16, st_sd == 1 ? Q16_ONE : 0, st_sd == 1 ? 0 : Q16_ONE, 0, false, dt_q16, rem);

      int32_t h_nv = q16_div(x_nv[1], x_nv[0] + x_nv[1]);
      if (st_nv == 1 && h_nv <= Q16_FROM_FLOAT(0.045) && x_nv[2] >= Q16_FROM_FLOAT(14.0) && x_nv[7] >= Q16_FROM_FLOAT(7.10)) {
        st_nv = 2;
      }
      rk4_cart8_step_q16_naive(x_nv, &cart_ws.params_q16, st_nv == 1 ? Q16_ONE : 0, st_nv == 1 ? 0 : Q16_ONE, 0, false, dt_q16);
    }

    double x_sd_f[8], x_nv_f[8];
    double norm_ref = 0.0, diff_sd = 0.0, diff_nv = 0.0, diff_32 = 0.0;
    for (int i = 0; i < 8; ++i) {
      x_sd_f[i] = Q16_TO_DOUBLE(x_sd[i]);
      x_nv_f[i] = Q16_TO_DOUBLE(x_nv[i]);
      double scale = (i == 2) ? 65.0 : ((i == 3) ? 25.0 : ((i == 7) ? 7.4 : 1.0));
      double dr = x64[i] / scale;
      double d_sd = (x_sd_f[i] - x64[i]) / scale;
      double d_nv = (x_nv_f[i] - x64[i]) / scale;
      double d_32 = ((double)x32[i] - x64[i]) / scale;
      norm_ref += dr * dr;
      diff_sd += d_sd * d_sd;
      diff_nv += d_nv * d_nv;
      diff_32 += d_32 * d_32;
    }
    err_sd_pct[s]    = 100.0 * sqrt(diff_sd / norm_ref);
    err_naive_pct[s] = 100.0 * sqrt(diff_nv / norm_ref);
    err_fp32_pct[s]  = 100.0 * sqrt(diff_32 / norm_ref);

    if (n_steps == 1000) {
      sw_fp64_1000   = sw64;
      sw_q16_1000    = sw_sd;
      n_fp64_1000    = x64[2];
      n_q16_1000     = x_sd_f[2];
      tscm_fp64_1000 = compute_tscm_pct(x64);
      tscm_q16_1000  = compute_tscm_pct(x_sd_f);
      h_fp64_1000    = 100.0 * x64[1] / (x64[0] + x64[1]);
      h_q16_1000     = 100.0 * x_sd_f[1] / (x_sd_f[0] + x_sd_f[1]);
    }
    assert(err_sd_pct[s] < 0.25);
  }

  /* ==========================================================================
   * Verify 100% Integer C99 Sequential 5-Channel Huber-EKF + NIS > 13.5
   * Biofouling Handover Under Donor Kinetic Mismatch & Electrode Fouling
   * ========================================================================== */
  memset(&cart_ws, 0, sizeof(cart_ws));
  init_cart_params_q16(&cart_ws.params_q16);

  /* Mismatched oncology donor plant (h0 = 49.0%, 0.90x mu_glyc, 1.15x k_Warburg, 0.88x k_mito_D, +4% OD bias) */
  double x_true[8] = {0.51, 0.49, 1.0, 24.5, 0.25, 1.20, 0.10, 7.32};
  const double x_nom0[8] = {0.54, 0.42, 1.0, 22.0, 0.28, 1.30, 0.08, 7.32};
  const double P0_diag[8] = {0.05, 0.05, 0.04, 4.0, 0.02, 0.08, 0.01, 0.01};
  const double Q_diag[8]  = {2.5e-4, 2.5e-4, 6e-3, 0.12, 2e-4, 8e-4, 1e-4, 2e-4};
  const double R_base_sig[5] = {0.55, 0.050, 0.35, 0.018, 0.012};
  double R_base_var[5];
  for (int m = 0; m < 5; ++m) {
    R_base_var[m] = R_base_sig[m] * R_base_sig[m];
    cart_ws.ekf_R_var_q20[m] = Q20_FROM_FLOAT(R_base_var[m]);
  }

  double x_hat64[8], P64[8][8];
  memset(P64, 0, sizeof(P64));
  for (int i = 0; i < 8; ++i) {
    x_hat64[i] = x_nom0[i];
    P64[i][i] = P0_diag[i];
    cart_ws.x_hat_q16[i] = Q16_FROM_FLOAT(x_nom0[i]);
    cart_ws.ekf_P_cov_q20[i][i] = Q20_FROM_FLOAT(P0_diag[i]);
    cart_ws.ekf_Q_diag_q20[i] = Q20_FROM_FLOAT(Q_diag[i]);
  }
  cart_ws.stage = 1U;
  int nis_streak64 = 0;
  bool ch0_dis64 = false;
  double handover_t_q16 = -1.0, handover_t_fp64 = -1.0;
  double gate_t_ekf_q16 = -1.0, gate_t_ekf_fp64 = -1.0;
  int quorum64 = 0;

  uint64_t rng_st = 20261111ULL;
  double ar_noise[5] = {0};
  const int n_ekf_steps = 560;
  const double dt_ekf = 168.0 / 560.0; /* 0.30 h = 18 min */
  const int32_t dt_ekf_q16 = Q16_FROM_FLOAT(dt_ekf);

  double sq_err_q16 = 0.0, sq_err_fp64 = 0.0, sq_diff_q16_vs_fp64 = 0.0;
  for (int k = 0; k < n_ekf_steps; ++k) {
    double t = k * dt_ekf;
    double u1 = (cart_ws.stage == 1U) ? 1.0 : 0.0;
    double u2 = (cart_ws.stage == 1U) ? 0.0 : 1.0;
    int32_t u1_q16 = (cart_ws.stage == 1U) ? Q16_ONE : 0;
    int32_t u2_q16 = (cart_ws.stage == 1U) ? 0 : Q16_ONE;

    double y_clean[5], y_meas[5];
    int32_t y_meas_q16[5];
    cassette_obs_fp64(x_true, u2, 0.365, 1.04, y_clean);
    for (int m = 0; m < 5; ++m) {
      double sig_m = R_base_sig[m];
      if (m == 1) sig_m = sqrt(0.050 * 0.050 + (0.045 * x_true[2]) * (0.045 * x_true[2]));
      if (m == 2) sig_m = sqrt(0.35 * 0.35 + (0.040 * x_true[2]) * (0.040 * x_true[2]));
      ar_noise[m] = 0.65 * ar_noise[m] + sqrt(1.0 - 0.65 * 0.65) * prng_normal01(&rng_st) * sig_m;
      y_meas[m] = y_clean[m] + ar_noise[m];
    }
    /* Inject progressive protein biofouling on Channel 0 (R_LP) starting at t = 38.0 h */
    if (t >= 38.0) {
      double spike = (prng_uniform01(&rng_st) < 0.14) ? ((prng_uniform01(&rng_st) < 0.5 ? -1.0 : 1.0) * 6.8) : 0.0;
      y_meas[0] += 0.42 * (t - 38.0) + spike;
      if (y_meas[0] < 1.5) y_meas[0] = 1.5;
    }
    for (int m = 0; m < 5; ++m) {
      y_meas_q16[m] = Q16_FROM_FLOAT(y_meas[m]);
    }

    bool newly_q16 = false, newly_fp64 = false;
    cart_ekf_step_q16(&cart_ws, y_meas_q16, u1_q16, u2_q16, dt_ekf_q16, &newly_q16);
    cart_ekf_step_fp64(x_hat64, P64, Q_diag, R_base_var, y_meas, u1, u2, dt_ekf, &nis_streak64, &ch0_dis64, &newly_fp64);
    if (newly_q16 && handover_t_q16 < 0.0) handover_t_q16 = t;
    if (newly_fp64 && handover_t_fp64 < 0.0) handover_t_fp64 = t;

    double h_true = x_true[1] / (x_true[0] + x_true[1]);
    double w_q16_d = Q16_TO_DOUBLE(cart_ws.x_hat_q16[0]);
    double d_q16_d = Q16_TO_DOUBLE(cart_ws.x_hat_q16[1]);
    double h_hat_q16_d = d_q16_d / ((w_q16_d + d_q16_d > 1e-6) ? (w_q16_d + d_q16_d) : 1e-6);
    double h_hat_fp64_d = x_hat64[1] / ((x_hat64[0] + x_hat64[1] > 1e-6) ? (x_hat64[0] + x_hat64[1]) : 1e-6);

    sq_err_q16 += (h_hat_q16_d - h_true) * (h_hat_q16_d - h_true);
    sq_err_fp64 += (h_hat_fp64_d - h_true) * (h_hat_fp64_d - h_true);
    sq_diff_q16_vs_fp64 += (h_hat_q16_d - h_hat_fp64_d) * (h_hat_q16_d - h_hat_fp64_d);

    /* Closed-loop persistence gating (4 frames = 1.2 h) */
    int32_t h_hat_q16 = q16_div(cart_ws.x_hat_q16[1], cart_ws.x_hat_q16[0] + cart_ws.x_hat_q16[1] + 64);
    bool cond_q16 = ((h_hat_q16 <= Q16_FROM_FLOAT(0.056) && cart_ws.x_hat_q16[2] >= Q16_FROM_FLOAT(18.0)) ||
                     (t >= 83.0 && h_hat_q16 <= Q16_FROM_FLOAT(0.072) && cart_ws.x_hat_q16[2] >= Q16_FROM_FLOAT(18.0))) &&
                    (y_meas_q16[3] >= Q16_FROM_FLOAT(7.10));
    cart_ws.quorum = (cart_ws.stage == 1U && cond_q16) ? (uint8_t)(cart_ws.quorum + 1U) : 0U;
    if (cart_ws.stage == 1U && cart_ws.quorum >= 4U) {
      cart_ws.stage = 2U;
      gate_t_ekf_q16 = t;
    }

    bool cond_fp64 = ((h_hat_fp64_d <= 0.056 && x_hat64[2] >= 18.0) ||
                      (t >= 83.0 && h_hat_fp64_d <= 0.072 && x_hat64[2] >= 18.0)) &&
                     (y_meas[3] >= 7.10);
    if (gate_t_ekf_fp64 < 0.0) {
      quorum64 = cond_fp64 ? (quorum64 + 1) : 0;
      if (quorum64 >= 4) gate_t_ekf_fp64 = t;
    }

    double u1_next = (cart_ws.stage == 1U) ? 1.0 : 0.0;
    double u2_next = (cart_ws.stage == 1U) ? 0.0 : 1.0;
    rk4_cart8_step_fp64_ex(x_true, u1_next, u2_next, 0.0, false, dt_ekf, 0.90, 0.88, 1.15);
  }

  double ekf_rmse_q16_pct  = 100.0 * sqrt(sq_err_q16 / (double)n_ekf_steps);
  double ekf_rmse_fp64_pct = 100.0 * sqrt(sq_err_fp64 / (double)n_ekf_steps);
  double ekf_q16_vs_fp64_rmse_pct = 100.0 * sqrt(sq_diff_q16_vs_fp64 / (double)n_ekf_steps);
  assert(cart_ws.ch0_disabled == 1U && ch0_dis64);
  assert(ekf_rmse_q16_pct < 3.5);
  assert(ekf_q16_vs_fp64_rmse_pct < 0.25);

  /* Verify Plant-Preserving Acidosis Flush in C99 Q16.16 (Hours 72-96 Occlusion) */
  int32_t xf_no[8], xf_rb[8];
  int64_t rem_no[8] = {0}, rem_rb[8] = {0};
  const double x0_init[8] = {0.54, 0.42, 1.0, 22.0, 0.28, 1.30, 0.08, 7.32};
  for (int i = 0; i < 8; ++i) {
    xf_no[i] = Q16_FROM_FLOAT(x0_init[i]);
    xf_rb[i] = xf_no[i];
  }
  int st_no = 1, st_rb = 1;
  double flush_until = -1.0;
  int32_t min_ph_no = Q16_FROM_FLOAT(9.0), min_ph_rb = Q16_FROM_FLOAT(9.0);
  int32_t dt1000_q16 = Q16_FROM_FLOAT(0.168);

  for (int k = 0; k < 1000; ++k) {
    double t = k * 0.168;
    bool fault = (t >= 72.0 && t <= 96.0);

    int32_t h_no = q16_div(xf_no[1], xf_no[0] + xf_no[1]);
    if (st_no == 1 && h_no <= Q16_FROM_FLOAT(0.045) && xf_no[2] >= Q16_FROM_FLOAT(14.0) && xf_no[7] >= Q16_FROM_FLOAT(7.10)) {
      st_no = 2;
    }
    rk4_cart8_step_q16_sd(xf_no, &cart_ws.params_q16, st_no == 1 ? Q16_ONE : 0, st_no == 1 ? 0 : Q16_ONE, 0, fault, dt1000_q16, rem_no);
    if (xf_no[7] < min_ph_no) min_ph_no = xf_no[7];

    int32_t h_rb = q16_div(xf_rb[1], xf_rb[0] + xf_rb[1]);
    if (st_rb == 1 && h_rb <= Q16_FROM_FLOAT(0.045) && xf_rb[2] >= Q16_FROM_FLOAT(14.0) && xf_rb[7] >= Q16_FROM_FLOAT(7.12) && t >= flush_until) {
      st_rb = 2;
    }
    if (st_rb == 2 && xf_rb[7] <= Q16_FROM_FLOAT(6.85)) {
      st_rb = 1;
      flush_until = t + 6.0;
    }
    int32_t u_fl = (t < flush_until || xf_rb[7] <= Q16_FROM_FLOAT(6.92)) ? Q16_ONE : 0;
    rk4_cart8_step_q16_sd(xf_rb, &cart_ws.params_q16, st_rb == 1 ? Q16_ONE : 0, st_rb == 1 ? 0 : Q16_ONE, u_fl, fault, dt1000_q16, rem_rb);
    if (xf_rb[7] < min_ph_rb) min_ph_rb = xf_rb[7];
  }

  double xf_no_d[8], xf_rb_d[8];
  for (int i = 0; i < 8; ++i) {
    xf_no_d[i] = Q16_TO_DOUBLE(xf_no[i]);
    xf_rb_d[i] = Q16_TO_DOUBLE(xf_rb[i]);
  }
  double tscm_fault_no = compute_tscm_pct(xf_no_d);
  double tscm_fault_rb = compute_tscm_pct(xf_rb_d);
  assert(Q16_TO_DOUBLE(min_ph_rb) > 6.84);
  assert(tscm_fault_rb > tscm_fault_no + 10.0);

  if (json_out) {
    printf("{\n");
    printf("  \"core_dtcm_workspace_bytes\": %zu,\n", core_bytes);
    printf("  \"cart_dtcm_workspace_bytes\": %zu,\n", cart_bytes);
    printf("  \"dtcm_bank0_capacity_bytes\": 8192,\n");
    printf("  \"cart_dtcm_bank0_utilization_pct\": %.2f,\n", 100.0 * (double)cart_bytes / 8192.0);
    printf("  \"n1000_switch_time_fp64_h\": %.3f,\n", sw_fp64_1000);
    printf("  \"n1000_switch_time_q16_h\": %.3f,\n", sw_q16_1000);
    printf("  \"n1000_expansion_fp64_fold\": %.4f,\n", n_fp64_1000);
    printf("  \"n1000_expansion_q16_fold\": %.4f,\n", n_q16_1000);
    printf("  \"n1000_tscm_fp64_pct\": %.4f,\n", tscm_fp64_1000);
    printf("  \"n1000_tscm_q16_pct\": %.4f,\n", tscm_q16_1000);
    printf("  \"n1000_h_fp64_pct\": %.4f,\n", h_fp64_1000);
    printf("  \"n1000_h_q16_pct\": %.4f,\n", h_q16_1000);
    printf("  \"step_sweep\": [\n");
    for (int s = 0; s < 5; ++s) {
      printf("    {\"steps\": %d, \"err_q16_sd_pct\": %.5f, \"err_q16_naive_pct\": %.5f, \"err_fp32_pct\": %.6f}%s\n",
             step_counts[s], err_sd_pct[s], err_naive_pct[s], err_fp32_pct[s], (s < 4) ? "," : "");
    }
    printf("  ],\n");
    printf("  \"c99_ekf_biofoul_verification\": {\n");
    printf("    \"ekf_h_rmse_fp64_pct\": %.4f,\n", ekf_rmse_fp64_pct);
    printf("    \"ekf_h_rmse_q16_pct\": %.4f,\n", ekf_rmse_q16_pct);
    printf("    \"ekf_q16_vs_fp64_rmse_pct\": %.5f,\n", ekf_q16_vs_fp64_rmse_pct);
    printf("    \"nis_handover_time_fp64_h\": %.2f,\n", handover_t_fp64);
    printf("    \"nis_handover_time_q16_h\": %.2f,\n", handover_t_q16);
    printf("    \"closed_loop_gate_fp64_h\": %.2f,\n", gate_t_ekf_fp64);
    printf("    \"closed_loop_gate_q16_h\": %.2f\n", gate_t_ekf_q16);
    printf("  },\n");
    printf("  \"c99_acidosis_fault\": {\n");
    printf("    \"unmitigated_min_ph\": %.4f,\n", Q16_TO_DOUBLE(min_ph_no));
    printf("    \"unmitigated_tscm_pct\": %.4f,\n", tscm_fault_no);
    printf("    \"bypass_flush_min_ph\": %.4f,\n", Q16_TO_DOUBLE(min_ph_rb));
    printf("    \"bypass_flush_tscm_pct\": %.4f\n", tscm_fault_rb);
    printf("  }\n");
    printf("}\n");
  } else {
    printf("[PASS] CAR-T Cassette Zero-Heap Static SRAM Workspace: %zu / 8192 B (%.2f%%)\n",
           cart_bytes, 100.0 * (double)cart_bytes / 8192.0);
    printf("[PASS] N=1000 Parity: sw=%.2fh (Q16=%.2fh), N=%.2fx (Q16=%.2fx), T_SCM=%.2f%% (Q16=%.2f%%), L2_err=%.4f%%\n",
           sw_fp64_1000, sw_q16_1000, n_fp64_1000, n_q16_1000, tscm_fp64_1000, tscm_q16_1000, err_sd_pct[1]);
    printf("[PASS] C99 Sequential 5-Ch Huber-EKF + NIS Handover: h_RMSE=%.3f%% (FP64=%.3f%%, diff=%.4f%%), Handover=%.1fh (FP64=%.1fh), Gate=%.1fh\n",
           ekf_rmse_q16_pct, ekf_rmse_fp64_pct, ekf_q16_vs_fp64_rmse_pct, handover_t_q16, handover_t_fp64, gate_t_ekf_q16);
    printf("[PASS] C99 Acidosis Bypass Flush: min_pH %.2f -> %.2f, T_SCM %.2f%% -> %.2f%%\n",
           Q16_TO_DOUBLE(min_ph_no), Q16_TO_DOUBLE(min_ph_rb), tscm_fault_no, tscm_fault_rb);
  }
  return 0;
}

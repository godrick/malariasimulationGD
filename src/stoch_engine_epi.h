// stoch_engine_epi.h
// Shared header for stochastic SEI tau-leaping engine.
// Defines StochEpiParams struct and step function declaration for use
// by both stoch_engine_epi.cpp and decoupled_loop_native.cpp.

#ifndef STOCH_ENGINE_EPI_H
#define STOCH_ENGINE_EPI_H

#include <vector>

struct StochEpiParams {
  // ---- Dimensions ----
  int nNodes, nM, nG, nE, nL, nP, nEIP, nStages, nPair, nState;
  int model_type;  // 0=Imperial decoupled, 1=SIS coupled, 2=SEIR coupled

  // ---- Index arrays (0-based) ----
  std::vector<int> egg_ix;   // nE * nG * nM
  std::vector<int> larv_ix;  // nL * nG * nM
  std::vector<int> pup_ix;   // nP * nG * nM
  std::vector<int> male_ix;  // nG * nM
  std::vector<int> unm_ix;   // nG * nM
  std::vector<int> fem_ix;   // nG * nG * nStages * nM  (4D flattened)

  // Human state indices (coupled SIS/SEIR, 0-based, -1 if N/A)
  std::vector<int> hS_ix, hI_ix, hE_ix, hR_ix;  // nNodes each
  std::vector<int> mosy_nodes;  // 0-based mosquito node indices (length nM)
  std::vector<int> human_nodes; // 0-based human node indices (length nH)
  int nH;

  // ---- Rates ----
  double rE, rL, rP, rEIP;
  double beta, nu;
  std::vector<double> beta_vec;

  // Per-node mortality (length nM)
  std::vector<double> muE, muL, muP, muM, muF;

  // Density dependence
  bool log_dd;
  std::vector<double> K, gamma_dd;  // length nM

  // ---- Genetics (from cube) ----
  std::vector<double> omega_inv;  // nG
  std::vector<double> phi, xiF, xiM;  // nG
  std::vector<double> eta;    // nG * nG, row-major: eta[i*nG + j]
  std::vector<double> B_mat;  // (nG*nG) * nG, column-major

  // ---- Movement (mosquito, nM x nM) ----
  bool has_move;
  std::vector<double> move_rates;  // nM
  std::vector<std::vector<int>> move_dest;  // per origin: list of dest indices
  std::vector<double> move_probs_flat;  // nM * nM, row-major (raw, diag=0)

  // ---- Human movement (coupled, nH x nH) ----
  bool has_hmove;
  std::vector<double> h_move_rates;  // nH
  std::vector<std::vector<int>> h_move_dest;
  std::vector<double> h_move_probs_flat;  // nH * nH, row-major

  // ---- Trap mortality ----
  std::vector<double> muM_node_base, muF_node_base;  // nM

  double tol;

  // ---- Infection parameters (SIS/SEIR) ----
  std::vector<double> c_vec, b_vec;  // nG
  std::vector<double> a_vec;  // nM (per-node biting rate)
  double muH_param, r_param, delta_param;

  // ---- Infection parameters (Imperial) ----
  std::vector<double> cT_vec, cD_vec, cU_vec;  // nG
  std::vector<double> W_age;  // na
  std::vector<double> fd;     // age-dependent (length na)
  double d1, ID0, kd, gamma1_imp;

  // ---- Human state pointer (Imperial decoupled) ----
  std::vector<double> human_state_buf;  // owned copy of human state
  const double* human_state_ptr;  // flat vector, na * 10
  int human_state_len;
  int human_na;

  // ---- Seasonal ----
  bool has_seasonal;
  int tmax_seasonal;
  std::vector<double> beta_tbl, nu_tbl, qE_tbl, qL_tbl, qP_tbl, qEIP_tbl;
  std::vector<double> muE_tbl, muL_tbl, muP_tbl, muM_tbl, muF_tbl;
  std::vector<double> K_tbl, gamma_tbl;
  std::vector<double> a_tbl, muH_tbl, r_tbl, delta_tbl;
  // Per-node seasonal flags: true if table length > tmax+1 (i.e. nM * (tmax+1))
  bool muE_is_vec, muL_is_vec, muP_is_vec, muM_is_vec, muF_is_vec;
  bool K_is_vec, gamma_is_vec;
  bool a_is_vec;

  // ---- Preallocated working memory ----
  std::vector<double> w_E, w_L, w_P, w_M, w_U;
  std::vector<double> w_F_stack;   // nPair * nStages * nM
  std::vector<double> w_dE, w_dL, w_dP, w_dM, w_dU;
  std::vector<double> w_dF_stack;  // nPair * nStages * nM
  std::vector<double> w_births, w_L_in, w_P_in;
  std::vector<double> w_M_emerge, w_F_emerge;
  std::vector<double> w_L_tot, w_dd_rate;
  std::vector<double> w_F_sum;     // nPair * nM
  std::vector<double> w_F_mort;    // nPair * nStages * nM (mortality events)
  std::vector<double> w_mate_p, w_denom_v, w_prob_buf;
  std::vector<int> w_valid_v, w_mates_out;
  std::vector<double> w_lambda_g;  // nG (Imperial infection rates)

  // ---- Incidence accumulator (cumulative new human infections per node) ----
  bool track_incidence;
  std::vector<int> inc_ix;  // 0-based state indices, length nHumanNodes
  int nHumanNodes;
};

// Allocate working memory (call after dimensions are set)
void stoch_epi_alloc_work(StochEpiParams* p);

// Perform one tau-leap step (modifies state in-place, calls R::rpois/rmultinom)
void stoch_epi_do_step(StochEpiParams* p, double* state, double t, double dt);

#endif  // STOCH_ENGINE_EPI_H

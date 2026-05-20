// det_engine_epi.cpp
// Native C++ ODE RHS for MGDrive2 mosquito SEI epidemiological model
//
// Supports:
// - Imperial decoupled: mosquito-only ODE, infection from external human state
// - SIS/SEIR coupled: mosquito + human ODE in one system
//
// Pattern A: Register with deSolve via extern "C" + R_registerRoutines.
// deSolve calls derivs directly via function pointer.

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "engine_common.h"
#include "solver.h"
using namespace Rcpp;
using namespace arma;

struct EpiEngine {
  // Dimensions
  int nNodes;     // total landscape nodes
  int nM;         // number of mosquito nodes
  int nG;         // genotypes
  int nE, nL, nP; // Erlang stages
  int nEIP;       // EIP Erlang stages
  int nStages;    // S + E1..EnEIP + I = nEIP + 2
  int nPair;      // nG * nG
  int nState;     // total state dimension (mosquito part only for Imperial; mosq+human for coupled)

  // Model type: 0=Imperial, 1=SIS, 2=SEIR
  int model_type;

  // Index arrays (0-based)
  arma::ucube egg_ix;    // (nE, nG, nM)
  arma::ucube larv_ix;   // (nL, nG, nM)
  arma::ucube pup_ix;    // (nP, nG, nM)
  arma::umat male_ix;    // (nG, nM)
  arma::umat unm_ix;     // (nG, nM)
  // fem_ix is 4D: (nG, nG, nStages, nM) — stored as a flat vector, accessed via idx4
  std::vector<unsigned int> fem_ix;  // length nG*nG*nStages*nM

  // Human state indices (only for coupled SIS/SEIR, 0-based, -1 if not applicable)
  std::vector<int> hS_ix, hI_ix, hE_ix, hR_ix;  // length nNodes each
  std::vector<int> mosy_nodes; // 0-based indices of mosquito nodes

  // Rates
  double rE, rL, rP, rEIP;
  double beta, nu;
  arma::vec beta_vec;

  // Per-node mortality (length nM)
  arma::vec muE, muL, muP, muM, muF;

  // Density dependence
  bool log_dd;
  arma::vec K, gamma;

  // Cube genetics
  arma::vec omega, omega_inv;
  arma::vec phi_cube, xiF, xiM;
  arma::mat eta;
  arma::mat B_mat;

  // Movement (mosquito)
  bool has_move;
  arma::sp_mat move_probs;
  arma::vec move_rates, move_rowsum;

  // Human movement (coupled only)
  bool has_hmove;
  arma::sp_mat h_move_probs;
  arma::vec h_move_rates, h_move_rowsum;
  std::vector<int> human_nodes; // 0-based indices of human nodes
  int nH;

  // Trap mortality
  arma::vec muM_node_base, muF_node_base;

  double tol;
  bool clamp_nonneg;

  // Infection parameters
  // SIS/SEIR: c_vec (length nG), b_vec (length nG), a (per-node nM), muH, r, delta
  arma::vec c_vec, b_vec;
  arma::vec a_vec;  // nM (per-node biting rate)
  double muH_param, r_param, delta_param;

  // Imperial: cT_vec, cD_vec, cU_vec (length nG), W_age (length na), d1, fd (length na), ID0, kd, gamma1
  arma::vec cT_vec, cD_vec, cU_vec;
  arma::vec W_age;  // precomputed (foi_age * av0 * rel_foi) / omega
  arma::vec fd;     // age-dependent fold-change in detection probability (length na)
  double d1, ID0, kd, gamma1_imp;

  // Seasonal support
  bool has_seasonal;
  int tmax_seasonal;
  std::vector<double> beta_tbl, nu_tbl, qE_tbl, qL_tbl, qP_tbl, qEIP_tbl;
  std::vector<double> muE_tbl, muL_tbl, muP_tbl, muM_tbl, muF_tbl;
  std::vector<double> K_tbl_s, gamma_tbl_s;
  std::vector<double> a_tbl, muH_tbl, r_tbl_s, delta_tbl;
  bool muE_is_vec, muL_is_vec, muP_is_vec, muM_is_vec, muF_is_vec;
  bool K_is_vec_s, gamma_is_vec_s;
  bool a_is_vec;

  // Preallocated working memory
  mutable arma::mat E_w, L_w, P_w, M_w, U_w;
  mutable arma::mat dE_w, dL_w, dP_w, dM_w, dU_w;
  mutable arma::mat F_stack, dF_stack;  // (nPair*nStages, nM)
  mutable arma::mat F_sum;              // (nPair, nM)
  mutable arma::mat births_w;
  mutable arma::vec births_flat_w, L_tot_w, dd_rate_w;
  mutable arma::vec m_vec_w, p_vec_w, u_vec_w, female_total_w, female_unmated_w, mate_total_w;
  mutable arma::mat W_w, mate_p_w;
  mutable arma::vec denom_w;
  mutable arma::uvec valid_w;
  mutable arma::mat YM_w, YU_w, YF_w, inbound_M_w, inbound_U_w, inbound_F_w;
};

// ---- Helper: access fem_ix as 4D ----
static inline unsigned int epi_fem_ix(const EpiEngine* eng, int gf, int gm, int stage, int node) {
  return eng->fem_ix[gf + gm * eng->nG + stage * eng->nPair + node * eng->nPair * eng->nStages];
}

// ---- Create engine ----
// [[Rcpp::export]]
SEXP epi_engine_create(
    int nNodes, int nM, int nG, int nE, int nL, int nP, int nEIP,
    int model_type,
    arma::ucube egg_ix, arma::ucube larv_ix, arma::ucube pup_ix,
    arma::umat male_ix, arma::umat unm_ix,
    Rcpp::IntegerVector fem_ix_flat,
    double rE, double rL, double rP, double rEIP,
    arma::vec muE, arma::vec muL, arma::vec muP, arma::vec muM, arma::vec muF,
    bool log_dd, arma::vec K, arma::vec gamma_dd,
    double beta, double nu,
    arma::vec omega, arma::vec phi_cube, arma::vec xiF, arma::vec xiM,
    arma::mat eta, arma::mat B_mat,
    bool has_move, arma::mat move_probs_dense, arma::vec move_rates,
    arma::vec muM_node_base, arma::vec muF_node_base,
    double tol,
    Rcpp::Nullable<Rcpp::NumericVector> c_vec_r,
    Rcpp::Nullable<Rcpp::NumericVector> b_vec_r,
    arma::vec a_vec, double muH_param, double r_param, double delta_param,
    Rcpp::Nullable<Rcpp::NumericVector> cT_vec_r,
    Rcpp::Nullable<Rcpp::NumericVector> cD_vec_r,
    Rcpp::Nullable<Rcpp::NumericVector> cU_vec_r,
    Rcpp::Nullable<Rcpp::NumericVector> W_age_r,
    double d1, Rcpp::NumericVector fd_r, double ID0, double kd, double gamma1_imp,
    Rcpp::IntegerVector mosy_nodes_r,
    Rcpp::IntegerVector human_nodes_r,
    Rcpp::IntegerVector hS_ix_r,
    Rcpp::IntegerVector hI_ix_r,
    Rcpp::IntegerVector hE_ix_r,
    Rcpp::IntegerVector hR_ix_r,
    bool has_hmove,
    arma::mat h_move_probs_dense,
    arma::vec h_move_rates,
    int nState
) {
  EpiEngine* eng = new EpiEngine();

  eng->nNodes = nNodes;
  eng->nM = nM;
  eng->nG = nG;
  eng->nE = nE;
  eng->nL = nL;
  eng->nP = nP;
  eng->nEIP = nEIP;
  eng->nStages = nEIP + 2;
  eng->nPair = nG * nG;
  eng->nState = nState;
  eng->model_type = model_type;

  // Convert 1-based R indices to 0-based
  eng->egg_ix = egg_ix - 1;
  eng->larv_ix = larv_ix - 1;
  eng->pup_ix = pup_ix - 1;
  eng->male_ix = male_ix - 1;
  eng->unm_ix = unm_ix - 1;

  // fem_ix: already flat from R, convert to 0-based
  eng->fem_ix.resize(fem_ix_flat.size());
  for (int i = 0; i < fem_ix_flat.size(); i++) {
    eng->fem_ix[i] = static_cast<unsigned int>(fem_ix_flat[i] - 1);
  }

  eng->rE = rE;
  eng->rL = rL;
  eng->rP = rP;
  eng->rEIP = rEIP;
  eng->beta = beta;
  eng->beta_vec = arma::vec(nM, arma::fill::value(beta));
  eng->nu = nu;

  auto expand_vec = [&](arma::vec v, int n) -> arma::vec {
    if ((int)v.n_elem == 1) return arma::vec(n, arma::fill::value(v(0)));
    return v;
  };
  eng->muE = expand_vec(muE, nM);
  eng->muL = expand_vec(muL, nM);
  eng->muP = expand_vec(muP, nM);
  eng->muM = expand_vec(muM, nM);
  eng->muF = expand_vec(muF, nM);

  eng->log_dd = log_dd;
  eng->K = expand_vec(K, nM);
  eng->gamma = expand_vec(gamma_dd, nM);

  eng->omega = omega;
  eng->omega_inv.set_size(omega.n_elem);
  for (arma::uword i = 0; i < omega.n_elem; i++)
    eng->omega_inv(i) = (omega(i) == 0.0) ? 1e3 : 1.0 / omega(i);
  eng->phi_cube = phi_cube;
  eng->xiF = xiF;
  eng->xiM = xiM;
  eng->eta = eta;
  eng->B_mat = B_mat;

  eng->has_move = has_move;
  eng->move_rates = move_rates;
  if (has_move) {
    eng->move_rowsum = arma::sum(move_probs_dense, 1);
    eng->move_probs = arma::sp_mat(move_probs_dense);
  }

  eng->muM_node_base = muM_node_base;
  eng->muF_node_base = muF_node_base;
  eng->tol = tol;
  eng->clamp_nonneg = true;

  // Infection parameters
  eng->a_vec = a_vec;
  eng->muH_param = muH_param;
  eng->r_param = r_param;
  eng->delta_param = delta_param;

  if (c_vec_r.isNotNull()) eng->c_vec = Rcpp::as<arma::vec>(c_vec_r);
  if (b_vec_r.isNotNull()) eng->b_vec = Rcpp::as<arma::vec>(b_vec_r);
  if (cT_vec_r.isNotNull()) eng->cT_vec = Rcpp::as<arma::vec>(cT_vec_r);
  if (cD_vec_r.isNotNull()) eng->cD_vec = Rcpp::as<arma::vec>(cD_vec_r);
  if (cU_vec_r.isNotNull()) eng->cU_vec = Rcpp::as<arma::vec>(cU_vec_r);
  if (W_age_r.isNotNull()) eng->W_age = Rcpp::as<arma::vec>(W_age_r);

  eng->d1 = d1;
  eng->fd = Rcpp::as<arma::vec>(fd_r);
  eng->ID0 = ID0;
  eng->kd = kd;
  eng->gamma1_imp = gamma1_imp;

  // Node indices (convert to 0-based)
  eng->mosy_nodes.resize(mosy_nodes_r.size());
  for (int i = 0; i < mosy_nodes_r.size(); i++) eng->mosy_nodes[i] = mosy_nodes_r[i] - 1;

  eng->human_nodes.resize(human_nodes_r.size());
  for (int i = 0; i < human_nodes_r.size(); i++) eng->human_nodes[i] = human_nodes_r[i] - 1;
  eng->nH = (int)human_nodes_r.size();

  // Human state indices (convert to 0-based, -1 for NA)
  auto conv_hix = [](Rcpp::IntegerVector rv) -> std::vector<int> {
    std::vector<int> out(rv.size());
    for (int i = 0; i < rv.size(); i++) {
      out[i] = Rcpp::IntegerVector::is_na(rv[i]) ? -1 : rv[i] - 1;
    }
    return out;
  };
  eng->hS_ix = conv_hix(hS_ix_r);
  eng->hI_ix = conv_hix(hI_ix_r);
  eng->hE_ix = conv_hix(hE_ix_r);
  eng->hR_ix = conv_hix(hR_ix_r);

  eng->has_hmove = has_hmove;
  if (has_hmove) {
    eng->h_move_rowsum = arma::sum(h_move_probs_dense, 1);
    eng->h_move_probs = arma::sp_mat(h_move_probs_dense);
    eng->h_move_rates = h_move_rates;
  }

  eng->has_seasonal = false;
  eng->tmax_seasonal = 0;

  // Preallocate working memory
  int nGnM = nG * nM;
  eng->E_w.set_size(nE, nGnM);
  eng->L_w.set_size(nL, nGnM);
  eng->P_w.set_size(nP, nGnM);
  eng->M_w.set_size(nG, nM);
  eng->U_w.set_size(nG, nM);
  eng->dE_w.set_size(nE, nGnM);
  eng->dL_w.set_size(nL, nGnM);
  eng->dP_w.set_size(nP, nGnM);
  eng->dM_w.set_size(nG, nM);
  eng->dU_w.set_size(nG, nM);

  int nFS = eng->nPair * eng->nStages;
  eng->F_stack.set_size(nFS, nM);
  eng->dF_stack.set_size(nFS, nM);
  eng->F_sum.set_size(eng->nPair, nM);

  eng->births_w.set_size(nG, nM);
  eng->births_flat_w.set_size(nGnM);
  eng->L_tot_w.set_size(nM);
  eng->dd_rate_w.set_size(nM);

  eng->W_w.set_size(nG, nG);
  eng->mate_p_w.set_size(nG, nG);
  eng->m_vec_w.set_size(nG);
  eng->p_vec_w.set_size(nG);
  eng->u_vec_w.set_size(nG);
  eng->female_total_w.set_size(nG);
  eng->female_unmated_w.set_size(nG);
  eng->mate_total_w.set_size(nG);
  eng->denom_w.set_size(nG);
  eng->valid_w.set_size(nG);

  if (has_move) {
    eng->YM_w.set_size(nG, nM);
    eng->YU_w.set_size(nG, nM);
    eng->YF_w.set_size(nFS, nM);
    eng->inbound_M_w.set_size(nG, nM);
    eng->inbound_U_w.set_size(nG, nM);
    eng->inbound_F_w.set_size(nFS, nM);
  }

  XPtr<EpiEngine> ptr(eng, true);
  return ptr;
}

// ---- Set seasonal tables ----
// [[Rcpp::export]]
void epi_engine_set_seasonal(
    SEXP engine_ptr, int tmax_seasonal,
    Rcpp::Nullable<Rcpp::NumericVector> beta_tbl_r,
    Rcpp::Nullable<Rcpp::NumericVector> nu_tbl_r,
    Rcpp::Nullable<Rcpp::NumericVector> qE_tbl_r,
    Rcpp::Nullable<Rcpp::NumericVector> qL_tbl_r,
    Rcpp::Nullable<Rcpp::NumericVector> qP_tbl_r,
    Rcpp::Nullable<Rcpp::NumericVector> qEIP_tbl_r,
    Rcpp::Nullable<Rcpp::NumericVector> muE_tbl_r,
    Rcpp::Nullable<Rcpp::NumericVector> muL_tbl_r,
    Rcpp::Nullable<Rcpp::NumericVector> muP_tbl_r,
    Rcpp::Nullable<Rcpp::NumericVector> muM_tbl_r,
    Rcpp::Nullable<Rcpp::NumericVector> muF_tbl_r,
    Rcpp::Nullable<Rcpp::NumericVector> K_tbl_r,
    Rcpp::Nullable<Rcpp::NumericVector> gamma_tbl_r,
    Rcpp::Nullable<Rcpp::NumericVector> a_tbl_r,
    Rcpp::Nullable<Rcpp::NumericVector> muH_tbl_r,
    Rcpp::Nullable<Rcpp::NumericVector> r_tbl_r,
    Rcpp::Nullable<Rcpp::NumericVector> delta_tbl_r
) {
  XPtr<EpiEngine> eng(engine_ptr);
  int tlen = tmax_seasonal + 1;

  eng->has_seasonal = true;
  eng->tmax_seasonal = tmax_seasonal;

  auto copy_scalar_tbl = [&](Rcpp::Nullable<Rcpp::NumericVector>& src, std::vector<double>& dst, double static_val) {
    if (src.isNotNull()) {
      Rcpp::NumericVector v(src);
      dst.assign(v.begin(), v.end());
    } else {
      dst.assign(tlen, static_val);
    }
  };

  auto copy_tbl = [&](Rcpp::Nullable<Rcpp::NumericVector>& src, std::vector<double>& dst, double static_val, bool& is_vec) {
    if (src.isNotNull()) {
      Rcpp::NumericVector v(src);
      dst.assign(v.begin(), v.end());
      is_vec = ((int)dst.size() > tlen);
    } else {
      dst.assign(tlen, static_val);
      is_vec = false;
    }
  };

  copy_scalar_tbl(beta_tbl_r, eng->beta_tbl, eng->beta);
  copy_scalar_tbl(nu_tbl_r, eng->nu_tbl, eng->nu);
  double qE_s = eng->rE / eng->nE;
  double qL_s = eng->rL / eng->nL;
  double qP_s = eng->rP / eng->nP;
  double qEIP_s = eng->rEIP / eng->nEIP;
  copy_scalar_tbl(qE_tbl_r, eng->qE_tbl, qE_s);
  copy_scalar_tbl(qL_tbl_r, eng->qL_tbl, qL_s);
  copy_scalar_tbl(qP_tbl_r, eng->qP_tbl, qP_s);
  copy_scalar_tbl(qEIP_tbl_r, eng->qEIP_tbl, qEIP_s);

  copy_tbl(muE_tbl_r, eng->muE_tbl, eng->muE(0), eng->muE_is_vec);
  copy_tbl(muL_tbl_r, eng->muL_tbl, eng->muL(0), eng->muL_is_vec);
  copy_tbl(muP_tbl_r, eng->muP_tbl, eng->muP(0), eng->muP_is_vec);
  copy_tbl(muM_tbl_r, eng->muM_tbl, eng->muM(0), eng->muM_is_vec);
  copy_tbl(muF_tbl_r, eng->muF_tbl, eng->muF(0), eng->muF_is_vec);
  copy_tbl(K_tbl_r, eng->K_tbl_s, eng->K(0), eng->K_is_vec_s);
  copy_tbl(gamma_tbl_r, eng->gamma_tbl_s, eng->gamma(0), eng->gamma_is_vec_s);

  copy_tbl(a_tbl_r, eng->a_tbl, eng->a_vec(0), eng->a_is_vec);
  copy_scalar_tbl(muH_tbl_r, eng->muH_tbl, eng->muH_param);
  copy_scalar_tbl(r_tbl_r, eng->r_tbl_s, eng->r_param);
  copy_scalar_tbl(delta_tbl_r, eng->delta_tbl, eng->delta_param);
}

// ---- Global pointer for deSolve native interface ----
static EpiEngine* g_epi_engine_ptr = nullptr;

// Global human state buffer for Imperial decoupled (owns the data to prevent
// dangling pointers when the R-side temporary is garbage-collected).
static std::vector<double> g_human_state_buf;
static const double* g_human_state = nullptr;
static int g_human_state_len = 0;
static int g_human_na = 0;  // number of age classes

// [[Rcpp::export]]
void epi_engine_set_native(SEXP engine_ptr) {
  XPtr<EpiEngine> eng(engine_ptr);
  g_epi_engine_ptr = eng.get();
}

// [[Rcpp::export]]
void epi_engine_clear_native() {
  g_epi_engine_ptr = nullptr;
  g_human_state_buf.clear();
  g_human_state = nullptr;
  g_human_state_len = 0;
  g_human_na = 0;
}

// [[Rcpp::export]]
void epi_engine_set_runtime(
    SEXP engine_ptr,
    arma::vec beta,
    arma::vec muM,
    arma::vec muF,
    arma::vec K,
    arma::vec gamma_dd,
    arma::vec a_vec
) {
  XPtr<EpiEngine> eng(engine_ptr);

  auto expand_vec = [](arma::vec v, int n, const char* name) -> arma::vec {
    if ((int)v.n_elem == 1) return arma::vec(n, arma::fill::value(v(0)));
    if ((int)v.n_elem != n) {
      stop("`%s` must have length 1 or nM.", name);
    }
    return v;
  };

  arma::vec beta_v = expand_vec(beta, eng->nM, "beta");
  eng->beta = beta_v(0);
  eng->beta_vec = beta_v;
  eng->muM = expand_vec(muM, eng->nM, "muM");
  eng->muF = expand_vec(muF, eng->nM, "muF");
  eng->K = expand_vec(K, eng->nM, "K");
  eng->gamma = expand_vec(gamma_dd, eng->nM, "gamma_dd");
  eng->a_vec = expand_vec(a_vec, eng->nM, "a_vec");
  eng->has_seasonal = false;
}

// [[Rcpp::export]]
void epi_engine_set_human_state(Rcpp::NumericVector h_state, int na) {
  g_human_state_buf.assign(h_state.begin(), h_state.end());
  g_human_state = g_human_state_buf.data();
  g_human_state_len = (int)g_human_state_buf.size();
  g_human_na = na;
}

// ---- Core derivative computation ----
// Implements the mosquito SEI ODE RHS, optionally with coupled human dynamics.
// This is the common code used by both the Rcpp wrapper and the native deSolve interface.
static void epi_compute_derivs(EpiEngine* eng, double t, double* y, double* ydot) {
  const int nM = eng->nM;
  const int nG = eng->nG;
  const int nE = eng->nE;
  const int nL = eng->nL;
  const int nP = eng->nP;
  const int nEIP = eng->nEIP;
  const int nStages = eng->nStages;
  const int nPair = eng->nPair;
  const int nState = eng->nState;
  const int nNodes = eng->nNodes;

  // Clamp
  if (eng->clamp_nonneg) {
    for (int i = 0; i < nState; i++) {
      if (y[i] < 0.0) y[i] = 0.0;
    }
  }

  // ---- Resolve seasonal parameters ----
  double beta_t = 0.0, nu_t, rE, rL, rP, rEIP_t;
  double muH_t, r_t, delta_t;
  const double* a_arr;
  const double* beta_arr = nullptr;

  if (eng->has_seasonal) {
    int ti = seasonal_tidx(t, eng->tmax_seasonal);
    beta_t = eng->beta_tbl[ti];
    nu_t = eng->nu_tbl[ti];
    rE = eng->qE_tbl[ti] * nE;
    rL = eng->qL_tbl[ti] * nL;
    rP = eng->qP_tbl[ti] * nP;
    rEIP_t = eng->qEIP_tbl[ti] * nEIP;
    muH_t = eng->muH_tbl[ti];
    r_t = eng->r_tbl_s[ti];
    delta_t = eng->delta_tbl[ti];

    // Update node-varying mortality and biting rate
    auto update_nv = [&](arma::vec& dst, const std::vector<double>& tbl, bool is_vec, int n) {
      if (is_vec) {
        for (int i = 0; i < n; i++) dst(i) = tbl[i + ti * n];
      } else {
        double v = tbl[ti];
        for (int i = 0; i < n; i++) dst(i) = v;
      }
    };
    update_nv(eng->muE, eng->muE_tbl, eng->muE_is_vec, nM);
    update_nv(eng->muL, eng->muL_tbl, eng->muL_is_vec, nM);
    update_nv(eng->muP, eng->muP_tbl, eng->muP_is_vec, nM);
    update_nv(eng->muM, eng->muM_tbl, eng->muM_is_vec, nM);
    update_nv(eng->muF, eng->muF_tbl, eng->muF_is_vec, nM);
    if (eng->log_dd) update_nv(eng->K, eng->K_tbl_s, eng->K_is_vec_s, nM);
    else update_nv(eng->gamma, eng->gamma_tbl_s, eng->gamma_is_vec_s, nM);
    update_nv(eng->a_vec, eng->a_tbl, eng->a_is_vec, nM);
    a_arr = eng->a_vec.memptr();
  } else {
    if ((int)eng->beta_vec.n_elem == nM) {
      beta_arr = eng->beta_vec.memptr();
    } else {
      beta_t = eng->beta;
    }
    nu_t = eng->nu;
    rE = eng->rE;
    rL = eng->rL;
    rP = eng->rP;
    rEIP_t = eng->rEIP;
    muH_t = eng->muH_param;
    r_t = eng->r_param;
    delta_t = eng->delta_param;
    a_arr = eng->a_vec.memptr();
  }

  // ---- Extract state ----
  for (int k = 0; k < nM; k++) {
    for (int g = 0; g < nG; g++) {
      for (int e = 0; e < nE; e++)
        eng->E_w(e, g + k * nG) = y[eng->egg_ix(e, g, k)];
      for (int l = 0; l < nL; l++)
        eng->L_w(l, g + k * nG) = y[eng->larv_ix(l, g, k)];
      for (int p = 0; p < nP; p++)
        eng->P_w(p, g + k * nG) = y[eng->pup_ix(p, g, k)];
      eng->M_w(g, k) = y[eng->male_ix(g, k)];
      eng->U_w(g, k) = y[eng->unm_ix(g, k)];
    }
    for (int s = 0; s < nStages; s++) {
      for (int gm = 0; gm < nG; gm++) {
        for (int gf = 0; gf < nG; gf++) {
          int row = gf + gm * nG + s * nPair;
          eng->F_stack(row, k) = y[epi_fem_ix(eng, gf, gm, s, k)];
        }
      }
    }
  }

  // ---- Births (from total females across all infection stages) ----
  eng->F_sum.zeros();
  for (int s = 0; s < nStages; s++) {
    for (int k = 0; k < nM; k++) {
      for (int r = 0; r < nPair; r++) {
        eng->F_sum(r, k) += eng->F_stack(r + s * nPair, k);
      }
    }
  }
  eng->births_w.zeros();
  for (int k = 0; k < nM; k++) {
    double beta_k = (beta_arr == nullptr) ? beta_t : beta_arr[k];
    eng->births_w.col(k) = beta_k * eng->B_mat.t() * eng->F_sum.col(k);
  }
  eng->births_flat_w = arma::vectorise(eng->births_w);

  // ---- Egg dynamics ----
  eng->dE_w.zeros();
  for (int k = 0; k < nM; k++) {
    double muE_k = eng->muE(k);
    for (int g = 0; g < nG; g++) {
      int col = g + k * nG;
      for (int e = 0; e < nE; e++)
        eng->dE_w(e, col) = -muE_k * eng->E_w(e, col);
    }
  }
  if (nE > 1) {
    for (int e = 0; e < nE - 1; e++) {
      for (int col = 0; col < nG * nM; col++) {
        double flux = rE * eng->E_w(e, col);
        eng->dE_w(e, col) -= flux;
        eng->dE_w(e + 1, col) += flux;
      }
    }
  }
  for (int col = 0; col < nG * nM; col++) {
    eng->dE_w(nE - 1, col) -= rE * eng->E_w(nE - 1, col);
    eng->dE_w(0, col) += eng->births_flat_w(col);
  }

  // ---- Larvae dynamics ----
  eng->dL_w.zeros();
  for (int col = 0; col < nG * nM; col++)
    eng->dL_w(0, col) = rE * eng->E_w(nE - 1, col);
  if (nL > 1) {
    for (int l = 0; l < nL - 1; l++) {
      for (int col = 0; col < nG * nM; col++) {
        double flux = rL * eng->L_w(l, col);
        eng->dL_w(l, col) -= flux;
        eng->dL_w(l + 1, col) += flux;
      }
    }
  }
  for (int col = 0; col < nG * nM; col++)
    eng->dL_w(nL - 1, col) -= rL * eng->L_w(nL - 1, col);

  // Density dependence
  eng->L_tot_w.zeros();
  for (int k = 0; k < nM; k++) {
    double total = 0.0;
    for (int g = 0; g < nG; g++) {
      int col = g + k * nG;
      for (int l = 0; l < nL; l++) total += eng->L_w(l, col);
    }
    eng->L_tot_w(k) = total;
  }
  if (eng->log_dd) {
    for (int k = 0; k < nM; k++)
      eng->dd_rate_w(k) = eng->muL(k) * (1.0 + eng->L_tot_w(k) / eng->K(k));
  } else {
    for (int k = 0; k < nM; k++)
      eng->dd_rate_w(k) = eng->muL(k) + eng->gamma(k) * eng->L_tot_w(k);
  }
  for (int k = 0; k < nM; k++) {
    double dd = eng->dd_rate_w(k);
    for (int g = 0; g < nG; g++) {
      int col = g + k * nG;
      for (int l = 0; l < nL; l++)
        eng->dL_w(l, col) -= dd * eng->L_w(l, col);
    }
  }

  // ---- Pupae dynamics ----
  eng->dP_w.zeros();
  for (int k = 0; k < nM; k++) {
    double muP_k = eng->muP(k);
    for (int g = 0; g < nG; g++) {
      int col = g + k * nG;
      for (int p = 0; p < nP; p++)
        eng->dP_w(p, col) = -muP_k * eng->P_w(p, col);
    }
  }
  for (int col = 0; col < nG * nM; col++)
    eng->dP_w(0, col) += rL * eng->L_w(nL - 1, col);
  if (nP > 1) {
    for (int p = 0; p < nP - 1; p++) {
      for (int col = 0; col < nG * nM; col++) {
        double flux = rP * eng->P_w(p, col);
        eng->dP_w(p, col) -= flux;
        eng->dP_w(p + 1, col) += flux;
      }
    }
  }

  // ---- Adult dynamics ----
  eng->dM_w.zeros();
  eng->dU_w.zeros();
  eng->dF_stack.zeros();

  // Male emergence
  for (int k = 0; k < nM; k++) {
    for (int g = 0; g < nG; g++) {
      int col = g + k * nG;
      double PnP = eng->P_w(nP - 1, col);
      double emerge = rP * (1.0 - eng->phi_cube(g)) * eng->xiM(g) * PnP;
      eng->dM_w(g, k) += emerge;
      eng->dP_w(nP - 1, col) -= emerge;
    }
  }

  // Female emergence and mating — emerged females enter S stage (stage 0)
  double tol = eng->tol;
  for (int k = 0; k < nM; k++) {
    for (int g = 0; g < nG; g++) {
      eng->m_vec_w(g) = eng->M_w(g, k);
      eng->p_vec_w(g) = eng->P_w(nP - 1, g + k * nG);
      eng->u_vec_w(g) = eng->U_w(g, k);
    }
    double tot_m = arma::accu(eng->m_vec_w);

    if (tot_m > tol) {
      for (int i = 0; i < nG; i++)
        for (int j = 0; j < nG; j++)
          eng->W_w(i, j) = eng->eta(i, j) * eng->m_vec_w(j);
      for (int i = 0; i < nG; i++) {
        eng->denom_w(i) = arma::accu(eng->W_w.row(i));
        eng->valid_w(i) = (eng->denom_w(i) > tol) ? 1 : 0;
      }
      for (int i = 0; i < nG; i++) {
        if (eng->valid_w(i)) {
          for (int j = 0; j < nG; j++)
            eng->mate_p_w(i, j) = eng->W_w(i, j) / eng->denom_w(i);
        } else {
          for (int j = 0; j < nG; j++)
            eng->mate_p_w(i, j) = 0.0;
        }
      }

      // Female emergence into S stage
      for (int g = 0; g < nG; g++)
        eng->female_total_w(g) = eng->valid_w(g) ? rP * eng->phi_cube(g) * eng->xiF(g) * eng->p_vec_w(g) : 0.0;

      if (arma::accu(eng->female_total_w) > 0) {
        for (int gf = 0; gf < nG; gf++) {
          for (int gm = 0; gm < nG; gm++) {
            int row = gf + gm * nG;  // S stage = stage 0
            eng->dF_stack(row, k) += eng->mate_p_w(gf, gm) * eng->female_total_w(gf);
          }
          eng->dP_w(nP - 1, gf + k * nG) -= eng->female_total_w(gf);
        }
      }

      // Unmated mating
      for (int g = 0; g < nG; g++)
        eng->mate_total_w(g) = eng->valid_w(g) ? nu_t * eng->u_vec_w(g) : 0.0;
      if (arma::accu(eng->mate_total_w) > 0) {
        for (int gf = 0; gf < nG; gf++) {
          eng->dU_w(gf, k) -= eng->mate_total_w(gf);
          for (int gm = 0; gm < nG; gm++) {
            int row = gf + gm * nG;
            eng->dF_stack(row, k) += eng->mate_p_w(gf, gm) * eng->mate_total_w(gf);
          }
        }
      }
    } else {
      for (int g = 0; g < nG; g++)
        eng->female_unmated_w(g) = rP * eng->phi_cube(g) * eng->xiF(g) * eng->p_vec_w(g);
      if (arma::accu(eng->female_unmated_w) > 0) {
        for (int g = 0; g < nG; g++) {
          eng->dU_w(g, k) += eng->female_unmated_w(g);
          eng->dP_w(nP - 1, g + k * nG) -= eng->female_unmated_w(g);
        }
      }
    }
  }

  // ---- Movement (mosquito) ----
  if (eng->has_move) {
    for (int k = 0; k < nM; k++)
      for (int g = 0; g < nG; g++)
        eng->YM_w(g, k) = eng->M_w(g, k) * eng->move_rates(k);
    eng->inbound_M_w = eng->YM_w * eng->move_probs;
    for (int k = 0; k < nM; k++)
      for (int g = 0; g < nG; g++)
        eng->dM_w(g, k) += eng->inbound_M_w(g, k) - eng->YM_w(g, k) * eng->move_rowsum(k);

    for (int k = 0; k < nM; k++)
      for (int g = 0; g < nG; g++)
        eng->YU_w(g, k) = eng->U_w(g, k) * eng->move_rates(k);
    eng->inbound_U_w = eng->YU_w * eng->move_probs;
    for (int k = 0; k < nM; k++)
      for (int g = 0; g < nG; g++)
        eng->dU_w(g, k) += eng->inbound_U_w(g, k) - eng->YU_w(g, k) * eng->move_rowsum(k);

    // F_stack movement (all stages move together)
    int nFS = nPair * nStages;
    for (int k = 0; k < nM; k++)
      for (int r = 0; r < nFS; r++)
        eng->YF_w(r, k) = eng->F_stack(r, k) * eng->move_rates(k);
    eng->inbound_F_w = eng->YF_w * eng->move_probs;
    for (int k = 0; k < nM; k++)
      for (int r = 0; r < nFS; r++)
        eng->dF_stack(r, k) += eng->inbound_F_w(r, k) - eng->YF_w(r, k) * eng->move_rowsum(k);
  }

  // ---- Adult mortality (all SEI stages) ----
  for (int k = 0; k < nM; k++) {
    double muM_mod = eng->muM(k) * eng->muM_node_base(k);
    for (int g = 0; g < nG; g++)
      eng->dM_w(g, k) -= muM_mod * eng->M_w(g, k) * eng->omega_inv(g);
  }
  for (int k = 0; k < nM; k++) {
    double muF_mod = eng->muF(k) * eng->muF_node_base(k);
    for (int g = 0; g < nG; g++)
      eng->dU_w(g, k) -= muF_mod * eng->U_w(g, k) * eng->omega_inv(g);
  }
  // Female mortality across all SEI stages
  for (int k = 0; k < nM; k++) {
    double muF_mod = eng->muF(k) * eng->muF_node_base(k);
    for (int s = 0; s < nStages; s++) {
      for (int gf = 0; gf < nG; gf++) {
        for (int gm = 0; gm < nG; gm++) {
          int row = gf + gm * nG + s * nPair;
          eng->dF_stack(row, k) -= muF_mod * eng->F_stack(row, k) * eng->omega_inv(gf);
        }
      }
    }
  }

  // ---- Infection dynamics ----
  if (eng->model_type == 0) {
    // Imperial decoupled: infection rate from external human state via global pointer
    if (g_human_state != nullptr && g_human_na > 0) {
      int na = g_human_na;
      // Human state is a flat vector of length na * 10
      // Columns: S=0, T=1, D=2, A=3, U=4, P=5, ICA=6, IB=7, ID=8, IVA=9
      const double* hs = g_human_state;
      double d1 = eng->d1, ID0 = eng->ID0, kd = eng->kd, g1 = eng->gamma1_imp;
      int n_fd = (int)eng->fd.n_elem;

      double Sum_T = 0.0, Sum_D_star = 0.0, Sum_U_star = 0.0;
      for (int i = 0; i < na; i++) {
        double T_h = hs[i + na * 1];  // T column
        double D_h = hs[i + na * 2];  // D column
        double A_h = hs[i + na * 3];  // A column
        double U_h = hs[i + na * 4];  // U column
        double ID_h = hs[i + na * 8]; // ID column

        double fd_i = eng->fd(i % n_fd);
        double p_det = d1 + (1.0 - d1) / (1.0 + fd_i * std::pow(ID_h / ID0, kd));
        double p_term = std::pow(p_det, g1);
        double w = eng->W_age(i);

        Sum_T += T_h * w;
        Sum_D_star += (D_h + A_h * p_term) * w;
        Sum_U_star += (U_h + A_h * (1.0 - p_term)) * w;
      }

      // lambda_g per female genotype
      for (int gf = 0; gf < nG; gf++) {
        double lambda_gf = eng->cT_vec(gf) * Sum_T + eng->cD_vec(gf) * Sum_D_star + eng->cU_vec(gf) * Sum_U_star;
        // Apply to all mate combinations: S -> E1
        for (int gm = 0; gm < nG; gm++) {
          int row_S = gf + gm * nG;
          int row_E1 = gf + gm * nG + 1 * nPair;
          for (int k = 0; k < nM; k++) {
            double F_S = eng->F_stack(row_S, k);
            double inf_flux = lambda_gf * F_S;
            eng->dF_stack(row_S, k) -= inf_flux;
            eng->dF_stack(row_E1, k) += inf_flux;
          }
        }
      }
    }
  } else {
    // SIS/SEIR coupled: infection from human state in SPN
    // Read human state from state vector
    for (int k = 0; k < nM; k++) {
      int node = eng->mosy_nodes[k];
      double H_S = (eng->hS_ix[node] >= 0) ? y[eng->hS_ix[node]] : 0.0;
      double H_I = (eng->hI_ix[node] >= 0) ? y[eng->hI_ix[node]] : 0.0;
      double NH;
      if (eng->model_type == 2) {
        double H_E = (eng->hE_ix[node] >= 0) ? y[eng->hE_ix[node]] : 0.0;
        double H_R = (eng->hR_ix[node] >= 0) ? y[eng->hR_ix[node]] : 0.0;
        NH = H_S + H_E + H_I + H_R;
      } else {
        NH = H_S + H_I;
      }
      double prevH = (NH > tol) ? H_I / NH : 0.0;

      // Mosquito infection: S -> E1
      for (int gf = 0; gf < nG; gf++) {
        double lambda = eng->c_vec(gf) * a_arr[k] * prevH;
        for (int gm = 0; gm < nG; gm++) {
          int row_S = gf + gm * nG;
          int row_E1 = gf + gm * nG + 1 * nPair;
          double F_S = eng->F_stack(row_S, k);
          double inf_flux = lambda * F_S;
          eng->dF_stack(row_S, k) -= inf_flux;
          eng->dF_stack(row_E1, k) += inf_flux;
        }
      }
    }
  }

  // ---- EIP progression: E1->E2->...->EnEIP->I ----
  for (int ei = 0; ei < nEIP; ei++) {
    int s_from = 1 + ei;  // E stages are 1..nEIP
    int s_to = (ei < nEIP - 1) ? s_from + 1 : nStages - 1;  // last E goes to I
    for (int k = 0; k < nM; k++) {
      for (int pair = 0; pair < nPair; pair++) {
        int row_from = pair + s_from * nPair;
        int row_to = pair + s_to * nPair;
        double flux = rEIP_t * eng->F_stack(row_from, k);
        eng->dF_stack(row_from, k) -= flux;
        eng->dF_stack(row_to, k) += flux;
      }
    }
  }

  // ---- Coupled human dynamics (SIS/SEIR only) ----
  // Initialize human derivatives to zero
  std::vector<double> dH_S(nNodes, 0.0), dH_I(nNodes, 0.0);
  std::vector<double> dH_E(nNodes, 0.0), dH_R(nNodes, 0.0);

  if (eng->model_type > 0) {
    // Read human state
    std::vector<double> H_S(nNodes, 0.0), H_I(nNodes, 0.0);
    std::vector<double> H_E(nNodes, 0.0), H_R(nNodes, 0.0);
    std::vector<double> NH(nNodes, 0.0);

    for (int n = 0; n < nNodes; n++) {
      if (eng->hS_ix[n] >= 0) H_S[n] = y[eng->hS_ix[n]];
      if (eng->hI_ix[n] >= 0) H_I[n] = y[eng->hI_ix[n]];
      if (eng->model_type == 2) {
        if (eng->hE_ix[n] >= 0) H_E[n] = y[eng->hE_ix[n]];
        if (eng->hR_ix[n] >= 0) H_R[n] = y[eng->hR_ix[n]];
        NH[n] = H_S[n] + H_E[n] + H_I[n] + H_R[n];
      } else {
        NH[n] = H_S[n] + H_I[n];
      }
    }

    // Human births and death
    for (int n = 0; n < nNodes; n++) {
      dH_S[n] += muH_t * NH[n];
      dH_S[n] -= muH_t * H_S[n];
      dH_I[n] -= muH_t * H_I[n];
      if (eng->model_type == 2) {
        dH_E[n] -= muH_t * H_E[n];
        dH_R[n] -= muH_t * H_R[n];
      }
    }

    // Recovery
    if (eng->model_type == 1) {
      // SIS
      for (int n = 0; n < nNodes; n++) {
        double rec = r_t * H_I[n];
        dH_I[n] -= rec;
        dH_S[n] += rec;
      }
    } else {
      // SEIR
      for (int n = 0; n < nNodes; n++) {
        double lat = delta_t * H_E[n];
        dH_E[n] -= lat;
        dH_I[n] += lat;
        double rec = r_t * H_I[n];
        dH_I[n] -= rec;
        dH_R[n] += rec;
      }
    }

    // Human infection: sum of b_vec * I_V_female per node
    for (int k = 0; k < nM; k++) {
      int node = eng->mosy_nodes[k];
      double sum_bI = 0.0;
      // I stage is last stage
      int s_I = nStages - 1;
      for (int gf = 0; gf < nG; gf++) {
        double Iv_f = 0.0;
        for (int gm = 0; gm < nG; gm++) {
          Iv_f += eng->F_stack(gf + gm * nG + s_I * nPair, k);
        }
        sum_bI += eng->b_vec(gf) * Iv_f;
      }

      if (NH[node] > tol) {
        double inf_h = a_arr[k] * (H_S[node] / NH[node]) * sum_bI;
        dH_S[node] -= inf_h;
        if (eng->model_type == 1) dH_I[node] += inf_h;
        else dH_E[node] += inf_h;
      }
    }

    // Human movement
    if (eng->has_hmove) {
      int nHn = eng->nH;
      // Simple approach: iterate human nodes
      if (eng->model_type == 1) {
        arma::mat Hmat(2, nHn);
        for (int i = 0; i < nHn; i++) {
          int n = eng->human_nodes[i];
          Hmat(0, i) = H_S[n];
          Hmat(1, i) = H_I[n];
        }
        arma::mat YH(2, nHn);
        for (int i = 0; i < nHn; i++) {
          YH(0, i) = Hmat(0, i) * eng->h_move_rates(i);
          YH(1, i) = Hmat(1, i) * eng->h_move_rates(i);
        }
        arma::mat inbound_H = YH * eng->h_move_probs;
        for (int i = 0; i < nHn; i++) {
          int n = eng->human_nodes[i];
          dH_S[n] += inbound_H(0, i) - YH(0, i) * eng->h_move_rowsum(i);
          dH_I[n] += inbound_H(1, i) - YH(1, i) * eng->h_move_rowsum(i);
        }
      } else {
        arma::mat Hmat(4, nHn);
        for (int i = 0; i < nHn; i++) {
          int n = eng->human_nodes[i];
          Hmat(0, i) = H_S[n];
          Hmat(1, i) = H_E[n];
          Hmat(2, i) = H_I[n];
          Hmat(3, i) = H_R[n];
        }
        arma::mat YH(4, nHn);
        for (int i = 0; i < nHn; i++)
          for (int r = 0; r < 4; r++)
            YH(r, i) = Hmat(r, i) * eng->h_move_rates(i);
        arma::mat inbound_H = YH * eng->h_move_probs;
        for (int i = 0; i < nHn; i++) {
          int n = eng->human_nodes[i];
          dH_S[n] += inbound_H(0, i) - YH(0, i) * eng->h_move_rowsum(i);
          dH_E[n] += inbound_H(1, i) - YH(1, i) * eng->h_move_rowsum(i);
          dH_I[n] += inbound_H(2, i) - YH(2, i) * eng->h_move_rowsum(i);
          dH_R[n] += inbound_H(3, i) - YH(3, i) * eng->h_move_rowsum(i);
        }
      }
    }
  }

  // ---- Pack derivatives into output ----
  for (int i = 0; i < nState; i++) ydot[i] = 0.0;

  for (int k = 0; k < nM; k++) {
    for (int g = 0; g < nG; g++) {
      for (int e = 0; e < nE; e++)
        ydot[eng->egg_ix(e, g, k)] = eng->dE_w(e, g + k * nG);
      for (int l = 0; l < nL; l++)
        ydot[eng->larv_ix(l, g, k)] = eng->dL_w(l, g + k * nG);
      for (int p = 0; p < nP; p++)
        ydot[eng->pup_ix(p, g, k)] = eng->dP_w(p, g + k * nG);
      ydot[eng->male_ix(g, k)] = eng->dM_w(g, k);
      ydot[eng->unm_ix(g, k)] = eng->dU_w(g, k);
    }
    for (int s = 0; s < nStages; s++) {
      for (int gm = 0; gm < nG; gm++) {
        for (int gf = 0; gf < nG; gf++) {
          int row = gf + gm * nG + s * nPair;
          ydot[epi_fem_ix(eng, gf, gm, s, k)] = eng->dF_stack(row, k);
        }
      }
    }
  }

  // Human derivatives (coupled only)
  if (eng->model_type > 0) {
    for (int n = 0; n < nNodes; n++) {
      if (eng->hS_ix[n] >= 0) ydot[eng->hS_ix[n]] = dH_S[n];
      if (eng->hI_ix[n] >= 0) ydot[eng->hI_ix[n]] = dH_I[n];
      if (eng->model_type == 2) {
        if (eng->hE_ix[n] >= 0) ydot[eng->hE_ix[n]] = dH_E[n];
        if (eng->hR_ix[n] >= 0) ydot[eng->hR_ix[n]] = dH_R[n];
      }
    }
  }
}

// ---- Rcpp wrapper ----
// [[Rcpp::export]]
arma::vec epi_engine_dxdt(SEXP engine_ptr, double t, arma::vec state) {
  XPtr<EpiEngine> eng(engine_ptr);
  arma::vec deriv(eng->nState, arma::fill::zeros);
  epi_compute_derivs(eng.get(), t, state.memptr(), deriv.memptr());
  return deriv;
}

// [[Rcpp::export]]
Rcpp::XPtr<Solver> create_epi_solver(
    SEXP engine_ptr,
    std::vector<double> init,
    double r_tol,
    double a_tol,
    size_t max_steps
) {
  auto eqs = [engine_ptr](const state_t& x, state_t& dxdt, double t) {
    arma::vec state_view(const_cast<double*>(x.data()), x.size(), false, true);
    arma::vec deriv = epi_engine_dxdt(engine_ptr, t, state_view);
    dxdt.assign(deriv.begin(), deriv.end());
  };

  return Rcpp::XPtr<Solver>(
      new Solver(init, eqs, r_tol, a_tol, max_steps),
      true
  );
}

// ---- Native deSolve interface ----
extern "C" void epi_derivs_native(
    int* neq, double* t, double* y, double* ydot,
    double* yout, int* ip
) {
  EpiEngine* eng = g_epi_engine_ptr;
  if (!eng) {
    for (int i = 0; i < *neq; i++) ydot[i] = 0.0;
    return;
  }
  epi_compute_derivs(eng, *t, y, ydot);
}

// Native routine registration is consolidated in register_routines.cpp

// det_engine_ento.cpp
// Entomology-only deterministic ODE engine for carrier frequency prediction.
//
// Stripped-down version of det_engine_epi.cpp with:
// - No infection stages (nStages=1: mated females only, no S/E/I distinction)
// - No human state variables
// - No infection/EIP/human dynamics in the RHS
//
// State vector: eggs + larvae + pupae + males + unmated + mated_females(G×G)
// For 58 nodes, G=10, nE=2, nL=3, nP=2: ~9,860 state variables
// (vs 22,736 for the full epi engine with infection stages)

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "engine_common.h"
#include "solver.h"
using namespace Rcpp;
using namespace arma;

struct EntoEngine {
  int nNodes;
  int nG;
  int nE, nL, nP;
  int nPair;       // nG * nG
  int nState;

  // Seasonal carrying-capacity support
  bool has_seasonal;
  int tmax_seasonal;
  std::vector<double> K_tbl_s;
  bool K_is_vec_s;

  // Index arrays (0-based)
  arma::ucube egg_ix;    // (nE, nG, nNodes)
  arma::ucube larv_ix;   // (nL, nG, nNodes)
  arma::ucube pup_ix;    // (nP, nG, nNodes)
  arma::umat male_ix;    // (nG, nNodes)
  arma::umat unm_ix;     // (nG, nNodes)
  arma::umat fem_ix;     // (nPair, nNodes) — single stage, no infection

  // Rates
  double rE, rL, rP;
  double beta, nu;
  arma::vec beta_vec;

  // Per-node mortality
  arma::vec muE, muL, muP, muM, muF;

  // Density dependence
  bool log_dd;
  arma::vec K, gamma;

  // Cube genetics
  arma::vec omega_inv;
  arma::vec phi_cube, xiF, xiM;
  arma::mat eta;
  arma::mat B_mat;

  // Movement
  bool has_move;
  arma::sp_mat move_probs;
  arma::vec move_rates, move_rowsum;

  // Trap mortality
  arma::vec muM_node_base, muF_node_base;

  double tol;

  // Preallocated working memory
  mutable arma::mat E_w, L_w, P_w, M_w, U_w;
  mutable arma::mat dE_w, dL_w, dP_w, dM_w, dU_w;
  mutable arma::mat F_w, dF_w;           // (nPair, nNodes) — single stage
  mutable arma::mat births_w;
  mutable arma::vec births_flat_w, L_tot_w, dd_rate_w;
  mutable arma::vec m_vec_w, p_vec_w, u_vec_w, female_total_w, female_unmated_w, mate_total_w;
  mutable arma::mat W_w, mate_p_w;
  mutable arma::vec denom_w;
  mutable arma::uvec valid_w;
  mutable arma::mat YM_w, YF_w, inbound_M_w, inbound_F_w;
};

// ---- Helper: access fem_ix ----
static inline unsigned int ento_fem_ix(const EntoEngine* eng, int gf, int gm, int node) {
  return eng->fem_ix(gf + gm * eng->nG, node);
}

// ---- Create engine ----
// [[Rcpp::export]]
SEXP ento_engine_create(
    int nNodes, int nG, int nE, int nL, int nP,
    arma::ucube egg_ix, arma::ucube larv_ix, arma::ucube pup_ix,
    arma::umat male_ix, arma::umat unm_ix,
    arma::umat fem_ix,
    double rE, double rL, double rP,
    arma::vec muE, arma::vec muL, arma::vec muP, arma::vec muM, arma::vec muF,
    bool log_dd, arma::vec K, arma::vec gamma_dd,
    double beta, double nu,
    arma::vec omega, arma::vec phi_cube, arma::vec xiF, arma::vec xiM,
    arma::mat eta, arma::mat B_mat,
    bool has_move, arma::mat move_probs_dense, arma::vec move_rates,
    arma::vec muM_node_base, arma::vec muF_node_base,
    double tol,
    int nState
) {
  EntoEngine* eng = new EntoEngine();

  eng->nNodes = nNodes;
  eng->nG = nG;
  eng->nE = nE;
  eng->nL = nL;
  eng->nP = nP;
  eng->nPair = nG * nG;
  eng->nState = nState;
  eng->has_seasonal = false;
  eng->tmax_seasonal = 0;
  eng->K_is_vec_s = false;

  // Convert 1-based R indices to 0-based
  eng->egg_ix = egg_ix - 1;
  eng->larv_ix = larv_ix - 1;
  eng->pup_ix = pup_ix - 1;
  eng->male_ix = male_ix - 1;
  eng->unm_ix = unm_ix - 1;
  eng->fem_ix = fem_ix - 1;

  eng->rE = rE;
  eng->rL = rL;
  eng->rP = rP;
  eng->beta = beta;
  eng->beta_vec = arma::vec(nNodes, arma::fill::value(beta));
  eng->nu = nu;

  auto expand_vec = [&](arma::vec v, int n) -> arma::vec {
    if ((int)v.n_elem == 1) return arma::vec(n, arma::fill::value(v(0)));
    return v;
  };
  eng->muE = expand_vec(muE, nNodes);
  eng->muL = expand_vec(muL, nNodes);
  eng->muP = expand_vec(muP, nNodes);
  eng->muM = expand_vec(muM, nNodes);
  eng->muF = expand_vec(muF, nNodes);

  eng->log_dd = log_dd;
  eng->K = expand_vec(K, nNodes);
  eng->gamma = expand_vec(gamma_dd, nNodes);

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

  // Preallocate working memory
  int nGnN = nG * nNodes;
  eng->E_w.set_size(nE, nGnN);
  eng->L_w.set_size(nL, nGnN);
  eng->P_w.set_size(nP, nGnN);
  eng->M_w.set_size(nG, nNodes);
  eng->U_w.set_size(nG, nNodes);
  eng->dE_w.set_size(nE, nGnN);
  eng->dL_w.set_size(nL, nGnN);
  eng->dP_w.set_size(nP, nGnN);
  eng->dM_w.set_size(nG, nNodes);
  eng->dU_w.set_size(nG, nNodes);

  eng->F_w.set_size(eng->nPair, nNodes);
  eng->dF_w.set_size(eng->nPair, nNodes);

  eng->births_w.set_size(nG, nNodes);
  eng->births_flat_w.set_size(nGnN);
  eng->L_tot_w.set_size(nNodes);
  eng->dd_rate_w.set_size(nNodes);

  eng->m_vec_w.set_size(nG);
  eng->p_vec_w.set_size(nG);
  eng->u_vec_w.set_size(nG);
  eng->female_total_w.set_size(nG);
  eng->female_unmated_w.set_size(nG);
  eng->mate_total_w.set_size(nG);
  eng->W_w.set_size(nG, nG);
  eng->mate_p_w.set_size(nG, nG);
  eng->denom_w.set_size(nG);
  eng->valid_w.set_size(nG);

  eng->YM_w.set_size(nG, nNodes);
  eng->inbound_M_w.set_size(nG, nNodes);
  eng->YF_w.set_size(eng->nPair, nNodes);
  eng->inbound_F_w.set_size(eng->nPair, nNodes);

  XPtr<EntoEngine> ptr(eng, true);
  return ptr;
}

// ---- Set runtime parameters (called each timestep for K updates) ----
// [[Rcpp::export]]
void ento_engine_set_runtime(
    SEXP engine_ptr,
    arma::vec beta_vec,
    arma::vec muM,
    arma::vec muF,
    arma::vec K,
    arma::vec gamma_dd
) {
  XPtr<EntoEngine> eng(engine_ptr);
  eng->beta_vec = beta_vec;
  eng->muM = muM;
  eng->muF = muF;
  eng->K = K;
  eng->gamma = gamma_dd;
}

// ---- Set seasonal carrying-capacity table ----
// [[Rcpp::export]]
void ento_engine_set_seasonal_K(
    SEXP engine_ptr,
    int tmax_seasonal,
    Rcpp::Nullable<Rcpp::NumericVector> K_tbl_r
) {
  XPtr<EntoEngine> eng(engine_ptr);
  eng->has_seasonal = false;
  eng->tmax_seasonal = 0;
  eng->K_tbl_s.clear();
  eng->K_is_vec_s = false;

  if (K_tbl_r.isNull()) {
    return;
  }

  Rcpp::NumericVector K_tbl(K_tbl_r);
  eng->K_tbl_s.assign(K_tbl.begin(), K_tbl.end());
  eng->has_seasonal = true;
  eng->tmax_seasonal = tmax_seasonal;
  eng->K_is_vec_s = ((int)eng->K_tbl_s.size() > (tmax_seasonal + 1));
}

// ---- Core derivative computation ----
static void ento_compute_derivs(EntoEngine* eng, double t, double* y, double* ydot) {
  const int nN = eng->nNodes;
  const int nG = eng->nG;
  const int nE = eng->nE;
  const int nL = eng->nL;
  const int nP = eng->nP;
  const int nPair = eng->nPair;
  const int nState = eng->nState;

  // Clamp negatives
  for (int i = 0; i < nState; i++)
    if (y[i] < 0.0) y[i] = 0.0;

  double rE = eng->rE, rL = eng->rL, rP = eng->rP;
  double nu_t = eng->nu;
  const double* beta_arr = eng->beta_vec.memptr();

  if (eng->has_seasonal) {
    int ti = seasonal_tidx(t, eng->tmax_seasonal);
    if (eng->K_is_vec_s) {
      for (int k = 0; k < nN; k++) {
        eng->K(k) = eng->K_tbl_s[k + ti * nN];
      }
    } else {
      double K_t = eng->K_tbl_s[ti];
      for (int k = 0; k < nN; k++) {
        eng->K(k) = K_t;
      }
    }
  }

  // ---- Extract state ----
  for (int k = 0; k < nN; k++) {
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
    for (int gm = 0; gm < nG; gm++)
      for (int gf = 0; gf < nG; gf++)
        eng->F_w(gf + gm * nG, k) = y[ento_fem_ix(eng, gf, gm, k)];
  }

  // ---- Births ----
  eng->births_w.zeros();
  for (int k = 0; k < nN; k++)
    eng->births_w.col(k) = beta_arr[k] * eng->B_mat.t() * eng->F_w.col(k);
  eng->births_flat_w = arma::vectorise(eng->births_w);

  // ---- Egg dynamics ----
  eng->dE_w.zeros();
  for (int k = 0; k < nN; k++) {
    double muE_k = eng->muE(k);
    for (int g = 0; g < nG; g++) {
      int col = g + k * nG;
      for (int e = 0; e < nE; e++)
        eng->dE_w(e, col) = -muE_k * eng->E_w(e, col);
    }
  }
  if (nE > 1) {
    for (int e = 0; e < nE - 1; e++)
      for (int col = 0; col < nG * nN; col++) {
        double flux = rE * eng->E_w(e, col);
        eng->dE_w(e, col) -= flux;
        eng->dE_w(e + 1, col) += flux;
      }
  }
  for (int col = 0; col < nG * nN; col++) {
    eng->dE_w(nE - 1, col) -= rE * eng->E_w(nE - 1, col);
    eng->dE_w(0, col) += eng->births_flat_w(col);
  }

  // ---- Larvae dynamics ----
  eng->dL_w.zeros();
  for (int col = 0; col < nG * nN; col++)
    eng->dL_w(0, col) = rE * eng->E_w(nE - 1, col);
  if (nL > 1) {
    for (int l = 0; l < nL - 1; l++)
      for (int col = 0; col < nG * nN; col++) {
        double flux = rL * eng->L_w(l, col);
        eng->dL_w(l, col) -= flux;
        eng->dL_w(l + 1, col) += flux;
      }
  }
  for (int col = 0; col < nG * nN; col++)
    eng->dL_w(nL - 1, col) -= rL * eng->L_w(nL - 1, col);

  // Density dependence
  eng->L_tot_w.zeros();
  for (int k = 0; k < nN; k++) {
    double total = 0.0;
    for (int g = 0; g < nG; g++) {
      int col = g + k * nG;
      for (int l = 0; l < nL; l++) total += eng->L_w(l, col);
    }
    eng->L_tot_w(k) = total;
  }
  if (eng->log_dd) {
    for (int k = 0; k < nN; k++)
      eng->dd_rate_w(k) = eng->muL(k) * (1.0 + eng->L_tot_w(k) / eng->K(k));
  } else {
    for (int k = 0; k < nN; k++)
      eng->dd_rate_w(k) = eng->muL(k) + eng->gamma(k) * eng->L_tot_w(k);
  }
  for (int k = 0; k < nN; k++) {
    double dd = eng->dd_rate_w(k);
    for (int g = 0; g < nG; g++) {
      int col = g + k * nG;
      for (int l = 0; l < nL; l++)
        eng->dL_w(l, col) -= dd * eng->L_w(l, col);
    }
  }

  // ---- Pupae dynamics ----
  eng->dP_w.zeros();
  for (int k = 0; k < nN; k++) {
    double muP_k = eng->muP(k);
    for (int g = 0; g < nG; g++) {
      int col = g + k * nG;
      for (int p = 0; p < nP; p++)
        eng->dP_w(p, col) = -muP_k * eng->P_w(p, col);
    }
  }
  for (int col = 0; col < nG * nN; col++)
    eng->dP_w(0, col) += rL * eng->L_w(nL - 1, col);
  if (nP > 1) {
    for (int p = 0; p < nP - 1; p++)
      for (int col = 0; col < nG * nN; col++) {
        double flux = rP * eng->P_w(p, col);
        eng->dP_w(p, col) -= flux;
        eng->dP_w(p + 1, col) += flux;
      }
  }

  // ---- Adult dynamics ----
  eng->dM_w.zeros();
  eng->dU_w.zeros();
  eng->dF_w.zeros();

  // Male emergence
  for (int k = 0; k < nN; k++)
    for (int g = 0; g < nG; g++) {
      int col = g + k * nG;
      double PnP = eng->P_w(nP - 1, col);
      double emerge = rP * (1.0 - eng->phi_cube(g)) * eng->xiM(g) * PnP;
      eng->dM_w(g, k) += emerge;
      eng->dP_w(nP - 1, col) -= emerge;
    }

  // Female emergence and mating
  double tol = eng->tol;
  for (int k = 0; k < nN; k++) {
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

      // Female emergence (directly into mated)
      for (int g = 0; g < nG; g++)
        eng->female_total_w(g) = eng->valid_w(g) ? rP * eng->phi_cube(g) * eng->xiF(g) * eng->p_vec_w(g) : 0.0;

      if (arma::accu(eng->female_total_w) > 0) {
        for (int gf = 0; gf < nG; gf++) {
          for (int gm = 0; gm < nG; gm++)
            eng->dF_w(gf + gm * nG, k) += eng->mate_p_w(gf, gm) * eng->female_total_w(gf);
          eng->dP_w(nP - 1, gf + k * nG) -= eng->female_total_w(gf);
        }
      }

      // Unmated mating
      for (int g = 0; g < nG; g++)
        eng->mate_total_w(g) = eng->valid_w(g) ? nu_t * eng->u_vec_w(g) : 0.0;
      if (arma::accu(eng->mate_total_w) > 0) {
        for (int gf = 0; gf < nG; gf++) {
          eng->dU_w(gf, k) -= eng->mate_total_w(gf);
          for (int gm = 0; gm < nG; gm++)
            eng->dF_w(gf + gm * nG, k) += eng->mate_p_w(gf, gm) * eng->mate_total_w(gf);
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

  // ---- Movement ----
  if (eng->has_move) {
    for (int k = 0; k < nN; k++)
      for (int g = 0; g < nG; g++)
        eng->YM_w(g, k) = eng->M_w(g, k) * eng->move_rates(k);
    eng->inbound_M_w = eng->YM_w * eng->move_probs;
    for (int k = 0; k < nN; k++)
      for (int g = 0; g < nG; g++)
        eng->dM_w(g, k) += eng->inbound_M_w(g, k) - eng->YM_w(g, k) * eng->move_rowsum(k);

    for (int k = 0; k < nN; k++)
      for (int r = 0; r < nPair; r++)
        eng->YF_w(r, k) = eng->F_w(r, k) * eng->move_rates(k);
    eng->inbound_F_w = eng->YF_w * eng->move_probs;
    for (int k = 0; k < nN; k++)
      for (int r = 0; r < nPair; r++)
        eng->dF_w(r, k) += eng->inbound_F_w(r, k) - eng->YF_w(r, k) * eng->move_rowsum(k);
  }

  // ---- Adult mortality ----
  for (int k = 0; k < nN; k++) {
    double muM_mod = eng->muM(k) * eng->muM_node_base(k);
    for (int g = 0; g < nG; g++)
      eng->dM_w(g, k) -= muM_mod * eng->M_w(g, k) * eng->omega_inv(g);
  }
  for (int k = 0; k < nN; k++) {
    double muF_mod = eng->muF(k) * eng->muF_node_base(k);
    for (int g = 0; g < nG; g++)
      eng->dU_w(g, k) -= muF_mod * eng->U_w(g, k) * eng->omega_inv(g);
  }
  for (int k = 0; k < nN; k++) {
    double muF_mod = eng->muF(k) * eng->muF_node_base(k);
    for (int gf = 0; gf < nG; gf++)
      for (int gm = 0; gm < nG; gm++)
        eng->dF_w(gf + gm * nG, k) -= muF_mod * eng->F_w(gf + gm * nG, k) * eng->omega_inv(gf);
  }

  // ---- Pack derivatives ----
  for (int i = 0; i < nState; i++) ydot[i] = 0.0;

  for (int k = 0; k < nN; k++) {
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
    for (int gm = 0; gm < nG; gm++)
      for (int gf = 0; gf < nG; gf++)
        ydot[ento_fem_ix(eng, gf, gm, k)] = eng->dF_w(gf + gm * nG, k);
  }
}

// ---- Rcpp wrapper for dxdt ----
// [[Rcpp::export]]
arma::vec ento_engine_dxdt(SEXP engine_ptr, double t, arma::vec state) {
  XPtr<EntoEngine> eng(engine_ptr);
  arma::vec deriv(eng->nState, arma::fill::zeros);
  ento_compute_derivs(eng.get(), t, state.memptr(), deriv.memptr());
  return deriv;
}

// ---- Create Boost Dopri5 solver (kept for compatibility) ----
// [[Rcpp::export]]
Rcpp::XPtr<Solver> create_ento_solver(
    SEXP engine_ptr,
    std::vector<double> init,
    double r_tol,
    double a_tol,
    size_t max_steps
) {
  auto eqs = [engine_ptr](const state_t& x, state_t& dxdt, double t) {
    Rcpp::XPtr<EntoEngine> eng(engine_ptr);
    ento_compute_derivs(eng.get(), t,
                        const_cast<double*>(x.data()),
                        dxdt.data());
  };
  Rcpp::XPtr<Solver> solver(new Solver(init, eqs, r_tol, a_tol, max_steps), true);
  return solver;
}

// =====================================================
// Native deSolve C interface
// =====================================================
// Direct C-level integration with deSolve's ODE solvers (LSODA, Adams).
// Eliminates R/C++ boundary overhead on each derivative call.
// Same pattern as customMGDrive2's lifecycle_derivs_native.

static EntoEngine* g_ento_engine_ptr = nullptr;

// [[Rcpp::export]]
void ento_engine_set_native(SEXP engine_ptr) {
  XPtr<EntoEngine> eng(engine_ptr);
  g_ento_engine_ptr = eng.get();
}

// [[Rcpp::export]]
void ento_engine_clear_native() {
  g_ento_engine_ptr = nullptr;
}

// Native derivative function for deSolve
// Signature: void(int* neq, double* t, double* y, double* ydot, double* yout, int* ip)
extern "C" void ento_derivs_native(
    int* neq, double* t, double* y, double* ydot,
    double* yout, int* ip
) {
  EntoEngine* eng = g_ento_engine_ptr;
  if (!eng) {
    for (int i = 0; i < *neq; i++) ydot[i] = 0.0;
    return;
  }
  ento_compute_derivs(eng, *t, y, ydot);
}

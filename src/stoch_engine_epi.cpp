// stoch_engine_epi.cpp
// Native C++ stochastic tau-leaping simulation for MGDrive2 mosquito SEI model.
//
// Supports:
// - Imperial decoupled: mosquito-only stochastic, infection from external human state
// - SIS/SEIR coupled: mosquito + stochastic human dynamics in one system
//
// Pattern B: Entire simulation loop in C. R calls .Call() once.

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rmath.h>
#include "engine_common.h"
#include "stoch_engine_epi.h"
using namespace Rcpp;

// ---- Index helpers ----
// fem_ix 4D: fem_ix[gf + gm*nG + s*nPair + k*nPair*nStages]
static inline int efi(int gf, int gm, int s, int k, int nG, int nPair, int nStages) {
  return gf + gm * nG + s * nPair + k * nPair * nStages;
}
// F_stack layout: row = gf + gm*nG + s*nPair, col = k
// Flat index: row + k * nFS where nFS = nPair * nStages
static inline int fsi(int row, int k, int nFS) {
  return row + k * nFS;
}
// Egg/larva/pupa 3D: ix3(stage, geno, node, nStages, nG)
static inline int ix3s(int stage, int geno, int node, int d1, int d2) {
  return stage + geno * d1 + node * d1 * d2;
}

// ---- Allocate working memory ----
void stoch_epi_alloc_work(StochEpiParams* p) {
  int nM = p->nM, nG = p->nG, nE = p->nE, nL = p->nL, nP = p->nP;
  int nPair = p->nPair, nStages = p->nStages;
  int nFS = nPair * nStages;

  p->w_E.resize(nE * nG * nM);
  p->w_L.resize(nL * nG * nM);
  p->w_P.resize(nP * nG * nM);
  p->w_M.resize(nG * nM);
  p->w_U.resize(nG * nM);
  p->w_F_stack.resize(nFS * nM);
  p->w_dE.resize(nE * nG * nM);
  p->w_dL.resize(nL * nG * nM);
  p->w_dP.resize(nP * nG * nM);
  p->w_dM.resize(nG * nM);
  p->w_dU.resize(nG * nM);
  p->w_dF_stack.resize(nFS * nM);
  p->w_births.resize(nG * nM);
  p->w_L_in.resize(nG * nM);
  p->w_P_in.resize(nG * nM);
  p->w_M_emerge.resize(nG * nM);
  p->w_F_emerge.resize(nG * nM);
  p->w_L_tot.resize(nM);
  p->w_dd_rate.resize(nM);
  p->w_F_sum.resize(nPair * nM);
  p->w_F_mort.resize(nFS * nM);
  p->w_mate_p.resize(nG * nG);
  p->w_denom_v.resize(nG);
  p->w_prob_buf.resize(nG);
  p->w_valid_v.resize(nG);
  p->w_mates_out.resize(nG);
  p->w_lambda_g.resize(nG);
}

// ---- Imperial infection rate ----
static void calc_lambda_g(StochEpiParams* p, double* lambda_g) {
  if (p->human_state_ptr == nullptr || p->human_na <= 0) {
    for (int g = 0; g < p->nG; g++) lambda_g[g] = 0.0;
    return;
  }
  int na = p->human_na;
  const double* hs = p->human_state_ptr;
  double d1v = p->d1, ID0v = p->ID0, kdv = p->kd, g1 = p->gamma1_imp;
  int n_fd = (int)p->fd.size();

  double Sum_T = 0.0, Sum_D_star = 0.0, Sum_U_star = 0.0;
  for (int i = 0; i < na; i++) {
    double T_h = hs[i + na * 1];
    double D_h = hs[i + na * 2];
    double A_h = hs[i + na * 3];
    double U_h = hs[i + na * 4];
    double ID_h = hs[i + na * 8];
    double fd_i = p->fd[i % n_fd];
    double p_det = d1v + (1.0 - d1v) / (1.0 + fd_i * std::pow(ID_h / ID0v, kdv));
    double p_term = std::pow(p_det, g1);
    double w = p->W_age[i];
    Sum_T += T_h * w;
    Sum_D_star += (D_h + A_h * p_term) * w;
    Sum_U_star += (U_h + A_h * (1.0 - p_term)) * w;
  }
  for (int g = 0; g < p->nG; g++) {
    lambda_g[g] = p->cT_vec[g] * Sum_T + p->cD_vec[g] * Sum_D_star + p->cU_vec[g] * Sum_U_star;
  }
}

// ---- One tau-leap step ----
void stoch_epi_do_step(StochEpiParams* p, double* state, double t, double dt) {
  const int nM = p->nM, nG = p->nG, nE = p->nE, nL = p->nL, nP = p->nP;
  const int nEIP = p->nEIP, nStages = p->nStages, nPair = p->nPair;
  const int nFS = nPair * nStages;
  const int nState = p->nState;
  const int nNodes = p->nNodes;
  const double tol = p->tol;

  // Aliases to working memory
  double* E = p->w_E.data();
  double* L = p->w_L.data();
  double* P = p->w_P.data();
  double* M = p->w_M.data();
  double* U = p->w_U.data();
  double* FS = p->w_F_stack.data();
  double* dE = p->w_dE.data();
  double* dL = p->w_dL.data();
  double* dP = p->w_dP.data();
  double* dM = p->w_dM.data();
  double* dU = p->w_dU.data();
  double* dFS = p->w_dF_stack.data();
  double* births = p->w_births.data();
  double* L_in = p->w_L_in.data();
  double* P_in = p->w_P_in.data();
  double* M_emerge = p->w_M_emerge.data();
  double* F_emerge = p->w_F_emerge.data();
  double* L_tot = p->w_L_tot.data();
  double* dd_rate = p->w_dd_rate.data();
  double* F_sum = p->w_F_sum.data();
  double* F_mort = p->w_F_mort.data();

  // ---- Clamp non-negative ----
  for (int i = 0; i < nState; i++)
    if (state[i] < 0.0) state[i] = 0.0;

  // ---- Resolve seasonal parameters ----
  double beta_t = 0.0, nu_t, rE_t, rL_t, rP_t, rEIP_t;
  double muH_t, r_t, delta_t;
  const double* a_arr;
  const double* beta_arr = nullptr;
  if (p->has_seasonal) {
    int ti = seasonal_tidx(t, p->tmax_seasonal);
    beta_t = p->beta_tbl[ti];
    nu_t = p->nu_tbl[ti];
    rE_t = p->qE_tbl[ti] * nE;
    rL_t = p->qL_tbl[ti] * nL;
    rP_t = p->qP_tbl[ti] * nP;
    rEIP_t = p->qEIP_tbl[ti] * nEIP;
    muH_t = p->muH_tbl[ti];
    r_t = p->r_tbl[ti];
    delta_t = p->delta_tbl[ti];
    // Override per-node mortality, DD, and biting rate with seasonal values.
    // Supports both scalar tables (same value all nodes) and per-node
    // tables (stride layout: tbl[node + ti * nM]).
    auto update_nv = [&](std::vector<double>& dst, const std::vector<double>& tbl, bool is_vec, int n) {
      if (is_vec) {
        for (int i = 0; i < n; i++) dst[i] = tbl[i + ti * n];
      } else {
        double v = tbl[ti];
        for (int i = 0; i < n; i++) dst[i] = v;
      }
    };
    update_nv(p->muE, p->muE_tbl, p->muE_is_vec, nM);
    update_nv(p->muL, p->muL_tbl, p->muL_is_vec, nM);
    update_nv(p->muP, p->muP_tbl, p->muP_is_vec, nM);
    update_nv(p->muM, p->muM_tbl, p->muM_is_vec, nM);
    update_nv(p->muF, p->muF_tbl, p->muF_is_vec, nM);
    update_nv(p->K, p->K_tbl, p->K_is_vec, nM);
    update_nv(p->gamma_dd, p->gamma_tbl, p->gamma_is_vec, nM);
    update_nv(p->a_vec, p->a_tbl, p->a_is_vec, nM);
    a_arr = p->a_vec.data();
  } else {
    if ((int)p->beta_vec.size() == nM) {
      beta_arr = p->beta_vec.data();
    } else {
      beta_t = p->beta;
    }
    nu_t = p->nu;
    rE_t = p->rE;
    rL_t = p->rL;
    rP_t = p->rP;
    rEIP_t = p->rEIP;
    muH_t = p->muH_param;
    r_t = p->r_param;
    delta_t = p->delta_param;
    a_arr = p->a_vec.data();
  }

  // ---- Extract state into working arrays ----
  for (int i = 0; i < nE * nG * nM; i++) E[i] = state[p->egg_ix[i]];
  for (int i = 0; i < nL * nG * nM; i++) L[i] = state[p->larv_ix[i]];
  for (int i = 0; i < nP * nG * nM; i++) P[i] = state[p->pup_ix[i]];
  for (int i = 0; i < nG * nM; i++) M[i] = state[p->male_ix[i]];
  for (int i = 0; i < nG * nM; i++) U[i] = state[p->unm_ix[i]];
  for (int k = 0; k < nM; k++)
    for (int s = 0; s < nStages; s++)
      for (int gm = 0; gm < nG; gm++)
        for (int gf = 0; gf < nG; gf++) {
          int row = gf + gm * nG + s * nPair;
          FS[fsi(row, k, nFS)] = state[p->fem_ix[efi(gf, gm, s, k, nG, nPair, nStages)]];
        }

  // ---- F_sum: total females across all stages (for births) ----
  for (int i = 0; i < nPair * nM; i++) F_sum[i] = 0.0;
  for (int k = 0; k < nM; k++)
    for (int s = 0; s < nStages; s++)
      for (int r = 0; r < nPair; r++)
        F_sum[r + k * nPair] += FS[fsi(r + s * nPair, k, nFS)];

  // ---- Births ----
  for (int k = 0; k < nM; k++) {
    for (int g = 0; g < nG; g++) {
      double s = 0.0;
      for (int pair = 0; pair < nPair; pair++)
        s += p->B_mat[pair + g * nPair] * F_sum[pair + k * nPair];
      double beta_k = (beta_arr == nullptr) ? beta_t : beta_arr[k];
      double lam = beta_k * s * dt;
      births[g + k * nG] = (lam > 0) ? R::rpois(lam) : 0.0;
    }
  }

  // ---- Egg dynamics ----
  std::fill(dE, dE + nE * nG * nM, 0.0);
  std::fill(L_in, L_in + nG * nM, 0.0);
  for (int k = 0; k < nM; k++) {
    double muE_k = p->muE[k];
    for (int g = 0; g < nG; g++) {
      for (int e = 0; e < nE; e++) {
        int idx = ix3s(e, g, k, nE, nG);
        double val = E[idx];
        double mort = (val > 0 && muE_k > 0) ? R::rpois(muE_k * val * dt) : 0.0;
        double adv = (val > 0 && rE_t > 0) ? R::rpois(rE_t * val * dt) : 0.0;
        dE[idx] -= mort + adv;
        if (e < nE - 1) dE[ix3s(e + 1, g, k, nE, nG)] += adv;
        else L_in[g + k * nG] = adv;
      }
      dE[ix3s(0, g, k, nE, nG)] += births[g + k * nG];
    }
  }

  // ---- Larvae dynamics ----
  std::fill(dL, dL + nL * nG * nM, 0.0);
  std::fill(P_in, P_in + nG * nM, 0.0);
  for (int k = 0; k < nM; k++) {
    double total = 0.0;
    for (int g = 0; g < nG; g++)
      for (int l = 0; l < nL; l++)
        total += L[ix3s(l, g, k, nL, nG)];
    L_tot[k] = total;
  }
  for (int k = 0; k < nM; k++) {
    if (p->log_dd)
      dd_rate[k] = p->muL[k] * (1.0 + L_tot[k] / p->K[k]);
    else
      dd_rate[k] = p->muL[k] + p->gamma_dd[k] * L_tot[k];
  }
  for (int k = 0; k < nM; k++) {
    for (int g = 0; g < nG; g++) {
      for (int l = 0; l < nL; l++) {
        int idx = ix3s(l, g, k, nL, nG);
        double val = L[idx];
        double mort = (val > 0 && dd_rate[k] > 0) ? R::rpois(dd_rate[k] * val * dt) : 0.0;
        double adv = (val > 0 && rL_t > 0) ? R::rpois(rL_t * val * dt) : 0.0;
        dL[idx] -= mort + adv;
        if (l < nL - 1) dL[ix3s(l + 1, g, k, nL, nG)] += adv;
        else P_in[g + k * nG] = adv;
      }
      dL[ix3s(0, g, k, nL, nG)] += L_in[g + k * nG];
    }
  }

  // ---- Pupae dynamics ----
  std::fill(dP, dP + nP * nG * nM, 0.0);
  std::fill(M_emerge, M_emerge + nG * nM, 0.0);
  std::fill(F_emerge, F_emerge + nG * nM, 0.0);
  for (int k = 0; k < nM; k++) {
    double muP_k = p->muP[k];
    for (int g = 0; g < nG; g++) {
      // Non-last stages
      for (int pp = 0; pp < nP - 1; pp++) {
        int idx = ix3s(pp, g, k, nP, nG);
        double val = P[idx];
        double mort = (val > 0 && muP_k > 0) ? R::rpois(muP_k * val * dt) : 0.0;
        double adv = (val > 0 && rP_t > 0) ? R::rpois(rP_t * val * dt) : 0.0;
        dP[idx] -= mort + adv;
        dP[ix3s(pp + 1, g, k, nP, nG)] += adv;
      }
      // Last pupal stage: competing hazards
      int idx_last = ix3s(nP - 1, g, k, nP, nG);
      double p_val = P[idx_last];
      double mort = (p_val > 0 && muP_k > 0) ? R::rpois(muP_k * p_val * dt) : 0.0;
      double m_pot = (p_val > 0 && rP_t > 0)
        ? R::rpois(rP_t * (1.0 - p->phi[g]) * p->xiM[g] * p_val * dt) : 0.0;

      // Check if males available for mating
      double tot_m = 0.0;
      for (int gg = 0; gg < nG; gg++) tot_m += M[gg + k * nG];
      double f_pot;
      if (tot_m > tol) {
        double W_row_sum = 0.0;
        for (int j = 0; j < nG; j++)
          W_row_sum += p->eta[g * nG + j] * M[j + k * nG];
        bool gvalid = (W_row_sum > tol);
        double f_rate = gvalid ? (rP_t * p->phi[g] * p->xiF[g] * p_val) : 0.0;
        f_pot = (f_rate > 0) ? R::rpois(f_rate * dt) : 0.0;
      } else {
        f_pot = (p_val > 0 && rP_t > 0)
          ? R::rpois(rP_t * p->phi[g] * p->xiF[g] * p_val * dt) : 0.0;
      }
      dP[idx_last] -= mort + m_pot + f_pot;
      M_emerge[g + k * nG] = m_pot;
      F_emerge[g + k * nG] = f_pot;
      // Input from previous stage
      dP[ix3s(0, g, k, nP, nG)] += P_in[g + k * nG];
    }
  }

  // ---- Male dynamics ----
  std::fill(dM, dM + nG * nM, 0.0);
  for (int k = 0; k < nM; k++) {
    double muM_mod = p->muM[k] * p->muM_node_base[k];
    for (int g = 0; g < nG; g++) {
      int idx = g + k * nG;
      double val = M[idx];
      double mort = (val > 0 && muM_mod > 0)
        ? R::rpois(muM_mod * val * p->omega_inv[g] * dt) : 0.0;
      dM[idx] = M_emerge[idx] - mort;
    }
  }

  // ---- Female dynamics (unmated + mated) ----
  std::fill(dU, dU + nG * nM, 0.0);
  std::fill(dFS, dFS + nFS * nM, 0.0);

  for (int k = 0; k < nM; k++) {
    double muF_mod = p->muF[k] * p->muF_node_base[k];

    // Unmated female mortality
    for (int g = 0; g < nG; g++) {
      int idx = g + k * nG;
      double val = U[idx];
      double mort = (val > 0 && muF_mod > 0)
        ? R::rpois(muF_mod * val * p->omega_inv[g] * dt) : 0.0;
      dU[idx] -= mort;
    }

    double tot_m = 0.0;
    for (int g = 0; g < nG; g++) tot_m += M[g + k * nG];

    if (tot_m > tol) {
      // Compute mating probabilities
      double* mate_p = p->w_mate_p.data();
      double* denom_v = p->w_denom_v.data();
      int* valid_v = p->w_valid_v.data();
      double* prob_buf = p->w_prob_buf.data();
      int* mates_out = p->w_mates_out.data();

      for (int i = 0; i < nG; i++) {
        denom_v[i] = 0.0;
        for (int j = 0; j < nG; j++) {
          double w = p->eta[i * nG + j] * M[j + k * nG];
          mate_p[i * nG + j] = w;
          denom_v[i] += w;
        }
        valid_v[i] = (denom_v[i] > tol) ? 1 : 0;
        if (valid_v[i]) {
          for (int j = 0; j < nG; j++)
            mate_p[i * nG + j] /= denom_v[i];
        } else {
          for (int j = 0; j < nG; j++)
            mate_p[i * nG + j] = 0.0;
        }
      }

      // Emerging females -> mated females (enter S stage = stage 0)
      for (int gf = 0; gf < nG; gf++) {
        int n_emerge = static_cast<int>(F_emerge[gf + k * nG]);
        if (n_emerge > 0) {
          for (int j = 0; j < nG; j++) prob_buf[j] = mate_p[gf * nG + j];
          R::rmultinom(n_emerge, prob_buf, nG, mates_out);
          for (int gm = 0; gm < nG; gm++) {
            int row_S = gf + gm * nG;  // stage 0
            dFS[fsi(row_S, k, nFS)] += mates_out[gm];
          }
        }
      }

      // Unmated female mating
      for (int gf = 0; gf < nG; gf++) {
        if (!valid_v[gf]) continue;
        double u_val = U[gf + k * nG];
        double u_rate = nu_t * u_val;
        int u_mate = (u_rate > 0) ? static_cast<int>(R::rpois(u_rate * dt)) : 0;
        if (u_mate > 0) {
          dU[gf + k * nG] -= u_mate;
          for (int j = 0; j < nG; j++) prob_buf[j] = mate_p[gf * nG + j];
          R::rmultinom(u_mate, prob_buf, nG, mates_out);
          for (int gm = 0; gm < nG; gm++) {
            int row_S = gf + gm * nG;
            dFS[fsi(row_S, k, nFS)] += mates_out[gm];
          }
        }
      }
    } else {
      // No males: emerging females go to unmated
      for (int g = 0; g < nG; g++)
        dU[g + k * nG] += F_emerge[g + k * nG];
    }
  }

  // ---- Female mortality (all SEI stages) ----
  for (int k = 0; k < nM; k++) {
    double muF_mod = p->muF[k] * p->muF_node_base[k];
    for (int s = 0; s < nStages; s++) {
      for (int gf = 0; gf < nG; gf++) {
        for (int gm = 0; gm < nG; gm++) {
          int row = gf + gm * nG + s * nPair;
          int fi = fsi(row, k, nFS);
          double val = FS[fi];
          double mort = (val > 0 && muF_mod > 0)
            ? R::rpois(muF_mod * val * p->omega_inv[gf] * dt) : 0.0;
          F_mort[fi] = mort;
        }
      }
    }
  }

  // ---- Infection: S -> E1 ----
  if (p->model_type == 0) {
    // Imperial decoupled: infection from external human state
    double* lambda_g = p->w_lambda_g.data();
    calc_lambda_g(p, lambda_g);
    for (int k = 0; k < nM; k++) {
      for (int gf = 0; gf < nG; gf++) {
        for (int gm = 0; gm < nG; gm++) {
          int row_S = gf + gm * nG;
          int row_E1 = gf + gm * nG + 1 * nPair;
          double F_S = FS[fsi(row_S, k, nFS)];
          double lam = lambda_g[gf] * F_S * dt;
          int inf = (lam > 0) ? static_cast<int>(R::rpois(lam)) : 0;
          dFS[fsi(row_S, k, nFS)] -= inf;
          dFS[fsi(row_E1, k, nFS)] += inf;
        }
      }
    }
  } else {
    // Coupled SIS/SEIR: infection from human state in SPN
    for (int k = 0; k < nM; k++) {
      int node = p->mosy_nodes[k];
      double H_S_val = (p->hS_ix[node] >= 0) ? state[p->hS_ix[node]] : 0.0;
      double H_I_val = (p->hI_ix[node] >= 0) ? state[p->hI_ix[node]] : 0.0;
      double NH;
      if (p->model_type == 2) {
        double H_E_val = (p->hE_ix[node] >= 0) ? state[p->hE_ix[node]] : 0.0;
        double H_R_val = (p->hR_ix[node] >= 0) ? state[p->hR_ix[node]] : 0.0;
        NH = H_S_val + H_E_val + H_I_val + H_R_val;
      } else {
        NH = H_S_val + H_I_val;
      }
      double prevH = (NH > tol) ? H_I_val / NH : 0.0;

      for (int gf = 0; gf < nG; gf++) {
        double lambda = p->c_vec[gf] * a_arr[k] * prevH;
        for (int gm = 0; gm < nG; gm++) {
          int row_S = gf + gm * nG;
          int row_E1 = gf + gm * nG + 1 * nPair;
          double F_S = FS[fsi(row_S, k, nFS)];
          double lam = lambda * F_S * dt;
          int inf = (lam > 0) ? static_cast<int>(R::rpois(lam)) : 0;
          dFS[fsi(row_S, k, nFS)] -= inf;
          dFS[fsi(row_E1, k, nFS)] += inf;
        }
      }
    }
  }

  // ---- EIP progression: E1->E2->...->EnEIP->I ----
  for (int ei = 0; ei < nEIP; ei++) {
    int s_from = 1 + ei;
    int s_to = (ei < nEIP - 1) ? s_from + 1 : nStages - 1;
    for (int k = 0; k < nM; k++) {
      for (int pair = 0; pair < nPair; pair++) {
        int row_from = pair + s_from * nPair;
        double val = FS[fsi(row_from, k, nFS)];
        double prog = (val > 0 && rEIP_t > 0) ? R::rpois(rEIP_t * val * dt) : 0.0;
        int row_to = pair + s_to * nPair;
        dFS[fsi(row_from, k, nFS)] -= prog;
        dFS[fsi(row_to, k, nFS)] += prog;
      }
    }
  }

  // ---- Subtract female mortality from dF_stack ----
  for (int i = 0; i < nFS * nM; i++)
    dFS[i] -= F_mort[i];

  // ---- Mosquito movement ----
  if (p->has_move) {
    for (int origin = 0; origin < nM; origin++) {
      double mr = p->move_rates[origin];
      if (mr <= 0 || p->move_dest[origin].empty()) continue;
      const auto& dests = p->move_dest[origin];
      int nDest = (int)dests.size();

      // Males
      for (int g = 0; g < nG; g++) {
        double m_val = M[g + origin * nG];
        if (m_val <= 0) continue;
        for (int d = 0; d < nDest; d++) {
          double p_raw = p->move_probs_flat[origin * nM + dests[d]];
          double lam = m_val * mr * dt * p_raw;
          if (lam > 0) {
            int flow = static_cast<int>(R::rpois(lam));
            if (flow > 0) {
              dM[g + origin * nG] -= flow;
              dM[g + dests[d] * nG] += flow;
            }
          }
        }
      }

      // Unmated females
      for (int g = 0; g < nG; g++) {
        double u_val = U[g + origin * nG];
        if (u_val <= 0) continue;
        for (int d = 0; d < nDest; d++) {
          double p_raw = p->move_probs_flat[origin * nM + dests[d]];
          double lam = u_val * mr * dt * p_raw;
          if (lam > 0) {
            int flow = static_cast<int>(R::rpois(lam));
            if (flow > 0) {
              dU[g + origin * nG] -= flow;
              dU[g + dests[d] * nG] += flow;
            }
          }
        }
      }

      // Females (all SEI stages)
      for (int row = 0; row < nFS; row++) {
        double f_val = FS[fsi(row, origin, nFS)];
        if (f_val <= 0) continue;
        for (int d = 0; d < nDest; d++) {
          double p_raw = p->move_probs_flat[origin * nM + dests[d]];
          double lam = f_val * mr * dt * p_raw;
          if (lam > 0) {
            int flow = static_cast<int>(R::rpois(lam));
            if (flow > 0) {
              dFS[fsi(row, origin, nFS)] -= flow;
              dFS[fsi(row, dests[d], nFS)] += flow;
            }
          }
        }
      }
    }
  }

  // ---- Coupled human dynamics (SIS/SEIR only) ----
  if (p->model_type > 0) {
    std::vector<double> H_S(nNodes, 0.0), H_I(nNodes, 0.0);
    std::vector<double> H_E(nNodes, 0.0), H_R(nNodes, 0.0);
    std::vector<double> NH(nNodes, 0.0);
    std::vector<double> dH_S(nNodes, 0.0), dH_I(nNodes, 0.0);
    std::vector<double> dH_E(nNodes, 0.0), dH_R(nNodes, 0.0);

    for (int n = 0; n < nNodes; n++) {
      if (p->hS_ix[n] >= 0) H_S[n] = state[p->hS_ix[n]];
      if (p->hI_ix[n] >= 0) H_I[n] = state[p->hI_ix[n]];
      if (p->model_type == 2) {
        if (p->hE_ix[n] >= 0) H_E[n] = state[p->hE_ix[n]];
        if (p->hR_ix[n] >= 0) H_R[n] = state[p->hR_ix[n]];
        NH[n] = H_S[n] + H_E[n] + H_I[n] + H_R[n];
      } else {
        NH[n] = H_S[n] + H_I[n];
      }
    }

    // Human births + death from S
    for (int n = 0; n < nNodes; n++) {
      if (p->hS_ix[n] < 0) continue;
      double b = R::rpois(muH_t * NH[n] * dt);
      dH_S[n] += b;
      double dS = (H_S[n] > 0 && muH_t > 0) ? R::rpois(muH_t * H_S[n] * dt) : 0.0;
      dH_S[n] -= dS;
      double dI = (H_I[n] > 0 && muH_t > 0) ? R::rpois(muH_t * H_I[n] * dt) : 0.0;
      dH_I[n] -= dI;
      if (p->model_type == 2) {
        double dE = (H_E[n] > 0 && muH_t > 0) ? R::rpois(muH_t * H_E[n] * dt) : 0.0;
        dH_E[n] -= dE;
        double dR = (H_R[n] > 0 && muH_t > 0) ? R::rpois(muH_t * H_R[n] * dt) : 0.0;
        dH_R[n] -= dR;
      }
    }

    // Human infection
    std::vector<int> new_inf_accum(nNodes, 0);
    for (int k = 0; k < nM; k++) {
      int node = p->mosy_nodes[k];
      if (p->hS_ix[node] < 0) continue;
      // Sum b_vec weighted infectious females
      double sum_bI = 0.0;
      int s_I = nStages - 1;
      for (int gf = 0; gf < nG; gf++) {
        double Iv_f = 0.0;
        for (int gm = 0; gm < nG; gm++)
          Iv_f += FS[fsi(gf + gm * nG + s_I * nPair, k, nFS)];
        sum_bI += p->b_vec[gf] * Iv_f;
      }
      if (NH[node] > tol) {
        double inf_rate = a_arr[k] * (H_S[node] / NH[node]) * sum_bI;
        int inf = (inf_rate > 0) ? static_cast<int>(R::rpois(inf_rate * dt)) : 0;
        new_inf_accum[node] += inf;
        dH_S[node] -= inf;
        if (p->model_type == 1) dH_I[node] += inf;
        else dH_E[node] += inf;
      }
    }

    // Recovery / latency
    if (p->model_type == 1) {
      // SIS: recovery I -> S
      for (int n = 0; n < nNodes; n++) {
        if (p->hI_ix[n] < 0) continue;
        int rec = (H_I[n] > 0 && r_t > 0) ? static_cast<int>(R::rpois(r_t * H_I[n] * dt)) : 0;
        dH_I[n] -= rec;
        dH_S[n] += rec;
      }
    } else {
      // SEIR: latency E->I, recovery I->R
      for (int n = 0; n < nNodes; n++) {
        if (p->hE_ix[n] < 0) continue;
        int lat = (H_E[n] > 0 && delta_t > 0) ? static_cast<int>(R::rpois(delta_t * H_E[n] * dt)) : 0;
        dH_E[n] -= lat;
        dH_I[n] += lat;
        int rec = (H_I[n] > 0 && r_t > 0) ? static_cast<int>(R::rpois(r_t * H_I[n] * dt)) : 0;
        dH_I[n] -= rec;
        dH_R[n] += rec;
      }
    }

    // Human movement
    if (p->has_hmove) {
      int nH = p->nH;
      for (int origin = 0; origin < nH; origin++) {
        double hr = p->h_move_rates[origin];
        if (hr <= 0 || p->h_move_dest[origin].empty()) continue;
        const auto& hdests = p->h_move_dest[origin];
        int nHDest = (int)hdests.size();
        int n_origin = p->human_nodes[origin];

        // S
        if (H_S[n_origin] > 0) {
          for (int d = 0; d < nHDest; d++) {
            double p_raw = p->h_move_probs_flat[origin * nH + hdests[d]];
            double lam = H_S[n_origin] * hr * dt * p_raw;
            if (lam > 0) {
              int flow = static_cast<int>(R::rpois(lam));
              if (flow > 0) {
                dH_S[n_origin] -= flow;
                dH_S[p->human_nodes[hdests[d]]] += flow;
              }
            }
          }
        }
        // I
        if (H_I[n_origin] > 0) {
          for (int d = 0; d < nHDest; d++) {
            double p_raw = p->h_move_probs_flat[origin * nH + hdests[d]];
            double lam = H_I[n_origin] * hr * dt * p_raw;
            if (lam > 0) {
              int flow = static_cast<int>(R::rpois(lam));
              if (flow > 0) {
                dH_I[n_origin] -= flow;
                dH_I[p->human_nodes[hdests[d]]] += flow;
              }
            }
          }
        }
        if (p->model_type == 2) {
          // E
          if (H_E[n_origin] > 0) {
            for (int d = 0; d < nHDest; d++) {
              double p_raw = p->h_move_probs_flat[origin * nH + hdests[d]];
              double lam = H_E[n_origin] * hr * dt * p_raw;
              if (lam > 0) {
                int flow = static_cast<int>(R::rpois(lam));
                if (flow > 0) {
                  dH_E[n_origin] -= flow;
                  dH_E[p->human_nodes[hdests[d]]] += flow;
                }
              }
            }
          }
          // R
          if (H_R[n_origin] > 0) {
            for (int d = 0; d < nHDest; d++) {
              double p_raw = p->h_move_probs_flat[origin * nH + hdests[d]];
              double lam = H_R[n_origin] * hr * dt * p_raw;
              if (lam > 0) {
                int flow = static_cast<int>(R::rpois(lam));
                if (flow > 0) {
                  dH_R[n_origin] -= flow;
                  dH_R[p->human_nodes[hdests[d]]] += flow;
                }
              }
            }
          }
        }
      }
    }

    // Apply human deltas to state
    for (int n = 0; n < nNodes; n++) {
      if (p->hS_ix[n] >= 0) state[p->hS_ix[n]] += dH_S[n];
      if (p->hI_ix[n] >= 0) state[p->hI_ix[n]] += dH_I[n];
      if (p->model_type == 2) {
        if (p->hE_ix[n] >= 0) state[p->hE_ix[n]] += dH_E[n];
        if (p->hR_ix[n] >= 0) state[p->hR_ix[n]] += dH_R[n];
      }
    }

    // Accumulate incidence into extended state slots
    if (p->track_incidence) {
      for (int n = 0; n < p->nHumanNodes; n++) {
        state[p->inc_ix[n]] += new_inf_accum[p->human_nodes[n]];
      }
    }
  }

  // ---- Apply mosquito deltas ----
  for (int i = 0; i < nE * nG * nM; i++) state[p->egg_ix[i]] += dE[i];
  for (int i = 0; i < nL * nG * nM; i++) state[p->larv_ix[i]] += dL[i];
  for (int i = 0; i < nP * nG * nM; i++) state[p->pup_ix[i]] += dP[i];
  for (int i = 0; i < nG * nM; i++) state[p->male_ix[i]] += dM[i];
  for (int i = 0; i < nG * nM; i++) state[p->unm_ix[i]] += dU[i];
  for (int k = 0; k < nM; k++)
    for (int s = 0; s < nStages; s++)
      for (int gm = 0; gm < nG; gm++)
        for (int gf = 0; gf < nG; gf++) {
          int row = gf + gm * nG + s * nPair;
          state[p->fem_ix[efi(gf, gm, s, k, nG, nPair, nStages)]] += dFS[fsi(row, k, nFS)];
        }

  // ---- Clamp non-negative ----
  for (int i = 0; i < nState; i++)
    if (state[i] < 0.0) state[i] = 0.0;
}

// ============================================================================
// ---- Rcpp exports ----
// ============================================================================

// [[Rcpp::export]]
SEXP stoch_epi_engine_create(
    int nNodes, int nM, int nG, int nE, int nL, int nP, int nEIP,
    int model_type,
    Rcpp::IntegerVector egg_ix_r, Rcpp::IntegerVector larv_ix_r, Rcpp::IntegerVector pup_ix_r,
    Rcpp::IntegerVector male_ix_r, Rcpp::IntegerVector unm_ix_r, Rcpp::IntegerVector fem_ix_r,
    double rE, double rL, double rP, double rEIP,
    Rcpp::NumericVector muE_r, Rcpp::NumericVector muL_r, Rcpp::NumericVector muP_r,
    Rcpp::NumericVector muM_r, Rcpp::NumericVector muF_r,
    bool log_dd, Rcpp::NumericVector K_r, Rcpp::NumericVector gamma_dd_r,
    double beta, double nu,
    Rcpp::NumericVector omega_inv_r, Rcpp::NumericVector phi_r,
    Rcpp::NumericVector xiF_r, Rcpp::NumericVector xiM_r,
    Rcpp::NumericMatrix eta_r, Rcpp::NumericMatrix B_mat_r,
    bool has_move, Rcpp::NumericMatrix move_probs_r, Rcpp::NumericVector move_rates_r,
    Rcpp::NumericVector muM_node_base_r, Rcpp::NumericVector muF_node_base_r,
    double tol,
    Rcpp::Nullable<Rcpp::NumericVector> c_vec_r,
    Rcpp::Nullable<Rcpp::NumericVector> b_vec_r,
    Rcpp::NumericVector a_vec_r, double muH_param, double r_param, double delta_param,
    Rcpp::Nullable<Rcpp::NumericVector> cT_vec_r,
    Rcpp::Nullable<Rcpp::NumericVector> cD_vec_r,
    Rcpp::Nullable<Rcpp::NumericVector> cU_vec_r,
    Rcpp::Nullable<Rcpp::NumericVector> W_age_r,
    double d1, Rcpp::NumericVector fd_r, double ID0, double kd, double gamma1_imp,
    Rcpp::IntegerVector mosy_nodes_r,
    Rcpp::IntegerVector human_nodes_r,
    Rcpp::IntegerVector hS_ix_r, Rcpp::IntegerVector hI_ix_r,
    Rcpp::IntegerVector hE_ix_r, Rcpp::IntegerVector hR_ix_r,
    bool has_hmove, Rcpp::NumericMatrix h_move_probs_r, Rcpp::NumericVector h_move_rates_r,
    int nState
) {
  StochEpiParams* p = new StochEpiParams();
  p->nNodes = nNodes; p->nM = nM; p->nG = nG;
  p->nE = nE; p->nL = nL; p->nP = nP; p->nEIP = nEIP;
  p->nStages = nEIP + 2;
  p->nPair = nG * nG;
  p->nState = nState;
  p->model_type = model_type;

  p->rE = rE; p->rL = rL; p->rP = rP; p->rEIP = rEIP;
  p->beta = beta; p->beta_vec.assign(nM, beta); p->nu = nu;
  p->log_dd = log_dd; p->tol = tol;
  p->has_move = has_move; p->has_hmove = has_hmove;

  // Convert 1-based indices to 0-based
  auto conv = [](Rcpp::IntegerVector v) -> std::vector<int> {
    std::vector<int> out(v.size());
    for (int i = 0; i < v.size(); i++) out[i] = v[i] - 1;
    return out;
  };
  p->egg_ix = conv(egg_ix_r);
  p->larv_ix = conv(larv_ix_r);
  p->pup_ix = conv(pup_ix_r);
  p->male_ix = conv(male_ix_r);
  p->unm_ix = conv(unm_ix_r);
  p->fem_ix = conv(fem_ix_r);

  // Human indices (convert, -1 for NA)
  auto conv_hix = [](Rcpp::IntegerVector rv) -> std::vector<int> {
    std::vector<int> out(rv.size());
    for (int i = 0; i < rv.size(); i++)
      out[i] = Rcpp::IntegerVector::is_na(rv[i]) ? -1 : rv[i] - 1;
    return out;
  };
  p->hS_ix = conv_hix(hS_ix_r);
  p->hI_ix = conv_hix(hI_ix_r);
  p->hE_ix = conv_hix(hE_ix_r);
  p->hR_ix = conv_hix(hR_ix_r);

  p->mosy_nodes.resize(mosy_nodes_r.size());
  for (int i = 0; i < mosy_nodes_r.size(); i++) p->mosy_nodes[i] = mosy_nodes_r[i] - 1;
  p->human_nodes.resize(human_nodes_r.size());
  for (int i = 0; i < human_nodes_r.size(); i++) p->human_nodes[i] = human_nodes_r[i] - 1;
  p->nH = (int)human_nodes_r.size();

  auto to_vec = [](Rcpp::NumericVector v) -> std::vector<double> {
    return std::vector<double>(v.begin(), v.end());
  };
  p->muE = to_vec(muE_r); p->muL = to_vec(muL_r); p->muP = to_vec(muP_r);
  p->muM = to_vec(muM_r); p->muF = to_vec(muF_r);
  p->K = to_vec(K_r); p->gamma_dd = to_vec(gamma_dd_r);
  p->omega_inv = to_vec(omega_inv_r);
  p->phi = to_vec(phi_r); p->xiF = to_vec(xiF_r); p->xiM = to_vec(xiM_r);

  // eta: row-major
  p->eta.resize(nG * nG);
  for (int i = 0; i < nG; i++)
    for (int j = 0; j < nG; j++)
      p->eta[i * nG + j] = eta_r(i, j);

  // B_mat: column-major
  p->B_mat.resize(B_mat_r.nrow() * B_mat_r.ncol());
  for (int c = 0; c < B_mat_r.ncol(); c++)
    for (int r = 0; r < B_mat_r.nrow(); r++)
      p->B_mat[r + c * B_mat_r.nrow()] = B_mat_r(r, c);

  p->muM_node_base = to_vec(muM_node_base_r);
  p->muF_node_base = to_vec(muF_node_base_r);
  p->move_rates = to_vec(move_rates_r);

  // Mosquito movement (nM x nM)
  if (has_move) {
    p->move_probs_flat.resize(nM * nM);
    for (int i = 0; i < nM; i++)
      for (int j = 0; j < nM; j++)
        p->move_probs_flat[i * nM + j] = move_probs_r(i, j);
    p->move_dest.resize(nM);
    for (int i = 0; i < nM; i++) {
      for (int j = 0; j < nM; j++) {
        if (p->move_probs_flat[i * nM + j] > 0 && p->move_rates[i] > 0)
          p->move_dest[i].push_back(j);
      }
    }
  }

  // Human movement (nH x nH)
  if (has_hmove) {
    int nH = p->nH;
    p->h_move_rates = to_vec(h_move_rates_r);
    p->h_move_probs_flat.resize(nH * nH);
    for (int i = 0; i < nH; i++)
      for (int j = 0; j < nH; j++)
        p->h_move_probs_flat[i * nH + j] = h_move_probs_r(i, j);
    p->h_move_dest.resize(nH);
    for (int i = 0; i < nH; i++) {
      for (int j = 0; j < nH; j++) {
        if (p->h_move_probs_flat[i * nH + j] > 0 && p->h_move_rates[i] > 0)
          p->h_move_dest[i].push_back(j);
      }
    }
  }

  // Infection params
  p->a_vec = to_vec(a_vec_r); p->muH_param = muH_param;
  p->r_param = r_param; p->delta_param = delta_param;
  if (c_vec_r.isNotNull()) p->c_vec = to_vec(Rcpp::NumericVector(c_vec_r));
  if (b_vec_r.isNotNull()) p->b_vec = to_vec(Rcpp::NumericVector(b_vec_r));
  if (cT_vec_r.isNotNull()) p->cT_vec = to_vec(Rcpp::NumericVector(cT_vec_r));
  if (cD_vec_r.isNotNull()) p->cD_vec = to_vec(Rcpp::NumericVector(cD_vec_r));
  if (cU_vec_r.isNotNull()) p->cU_vec = to_vec(Rcpp::NumericVector(cU_vec_r));
  if (W_age_r.isNotNull()) p->W_age = to_vec(Rcpp::NumericVector(W_age_r));
  p->d1 = d1; p->fd = std::vector<double>(fd_r.begin(), fd_r.end()); p->ID0 = ID0; p->kd = kd; p->gamma1_imp = gamma1_imp;

  p->human_state_ptr = nullptr;
  p->human_state_len = 0;
  p->human_na = 0;
  p->has_seasonal = false;
  p->tmax_seasonal = 0;

  // Incidence accumulator (off by default; set via stoch_epi_engine_set_incidence)
  p->track_incidence = false;
  p->nHumanNodes = 0;

  stoch_epi_alloc_work(p);

  Rcpp::XPtr<StochEpiParams> ptr(p, true);
  return ptr;
}

// [[Rcpp::export]]
void stoch_epi_engine_set_seasonal(
    SEXP engine_ptr, int tmax,
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
  Rcpp::XPtr<StochEpiParams> p(engine_ptr);
  int tlen = tmax + 1;
  p->has_seasonal = true;
  p->tmax_seasonal = tmax;

  // Scalar-only table copy (rates, etc.)
  auto copy_scalar_tbl = [&](Rcpp::Nullable<Rcpp::NumericVector>& src, std::vector<double>& dst, double def) {
    if (src.isNotNull()) {
      Rcpp::NumericVector v(src);
      dst.assign(v.begin(), v.end());
    } else {
      dst.assign(tlen, def);
    }
  };

  // Per-node table copy: detects vector tables (length > tlen means nM * tlen layout)
  auto copy_nv_tbl = [&](Rcpp::Nullable<Rcpp::NumericVector>& src, std::vector<double>& dst, double def, bool& is_vec) {
    if (src.isNotNull()) {
      Rcpp::NumericVector v(src);
      dst.assign(v.begin(), v.end());
      is_vec = ((int)dst.size() > tlen);
    } else {
      dst.assign(tlen, def);
      is_vec = false;
    }
  };

  double qE_s = p->rE / p->nE, qL_s = p->rL / p->nL, qP_s = p->rP / p->nP;
  double qEIP_s = p->rEIP / p->nEIP;
  copy_scalar_tbl(beta_tbl_r, p->beta_tbl, p->beta);
  copy_scalar_tbl(nu_tbl_r, p->nu_tbl, p->nu);
  copy_scalar_tbl(qE_tbl_r, p->qE_tbl, qE_s);
  copy_scalar_tbl(qL_tbl_r, p->qL_tbl, qL_s);
  copy_scalar_tbl(qP_tbl_r, p->qP_tbl, qP_s);
  copy_scalar_tbl(qEIP_tbl_r, p->qEIP_tbl, qEIP_s);
  copy_nv_tbl(muE_tbl_r, p->muE_tbl, p->muE[0], p->muE_is_vec);
  copy_nv_tbl(muL_tbl_r, p->muL_tbl, p->muL[0], p->muL_is_vec);
  copy_nv_tbl(muP_tbl_r, p->muP_tbl, p->muP[0], p->muP_is_vec);
  copy_nv_tbl(muM_tbl_r, p->muM_tbl, p->muM[0], p->muM_is_vec);
  copy_nv_tbl(muF_tbl_r, p->muF_tbl, p->muF[0], p->muF_is_vec);
  copy_nv_tbl(K_tbl_r, p->K_tbl, p->K[0], p->K_is_vec);
  copy_nv_tbl(gamma_tbl_r, p->gamma_tbl, p->gamma_dd[0], p->gamma_is_vec);
  copy_nv_tbl(a_tbl_r, p->a_tbl, p->a_vec[0], p->a_is_vec);
  copy_scalar_tbl(muH_tbl_r, p->muH_tbl, p->muH_param);
  copy_scalar_tbl(r_tbl_r, p->r_tbl, p->r_param);
  copy_scalar_tbl(delta_tbl_r, p->delta_tbl, p->delta_param);
}

// [[Rcpp::export]]
void stoch_epi_engine_set_runtime(
    SEXP engine_ptr,
    Rcpp::NumericVector beta_r,
    Rcpp::NumericVector muM_r,
    Rcpp::NumericVector muF_r,
    Rcpp::NumericVector K_r,
    Rcpp::NumericVector gamma_dd_r,
    Rcpp::NumericVector a_vec_r
) {
  Rcpp::XPtr<StochEpiParams> p(engine_ptr);

  auto expand_vec = [&](Rcpp::NumericVector src, int n, const char* name) -> std::vector<double> {
    if (src.size() == 1) return std::vector<double>(n, src[0]);
    if (src.size() != n) {
      stop("`%s` must have length 1 or nM.", name);
    }
    return std::vector<double>(src.begin(), src.end());
  };

  std::vector<double> beta_v = expand_vec(beta_r, p->nM, "beta");
  p->beta = beta_v[0];
  p->beta_vec = beta_v;
  p->muM = expand_vec(muM_r, p->nM, "muM");
  p->muF = expand_vec(muF_r, p->nM, "muF");
  p->K = expand_vec(K_r, p->nM, "K");
  p->gamma_dd = expand_vec(gamma_dd_r, p->nM, "gamma_dd");
  p->a_vec = expand_vec(a_vec_r, p->nM, "a_vec");
  p->has_seasonal = false;
}

// [[Rcpp::export]]
void stoch_epi_engine_set_incidence(SEXP engine_ptr, Rcpp::IntegerVector inc_ix_r) {
  Rcpp::XPtr<StochEpiParams> p(engine_ptr);
  p->nHumanNodes = (int)inc_ix_r.size();
  p->track_incidence = (p->nHumanNodes > 0);
  p->inc_ix.resize(p->nHumanNodes);
  for (int i = 0; i < p->nHumanNodes; i++) p->inc_ix[i] = inc_ix_r[i] - 1;
}

// [[Rcpp::export]]
void stoch_epi_set_human_state(SEXP engine_ptr, Rcpp::NumericVector h_state, int na) {
  Rcpp::XPtr<StochEpiParams> p(engine_ptr);
  p->human_state_buf.assign(h_state.begin(), h_state.end());
  p->human_state_ptr = p->human_state_buf.data();
  p->human_state_len = (int)p->human_state_buf.size();
  p->human_na = na;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix stoch_epi_simulate(
    SEXP engine_ptr,
    Rcpp::NumericVector x0,
    double dt_stoch,
    double dt_out,
    double tmax
) {
  Rcpp::XPtr<StochEpiParams> p(engine_ptr);
  const int nState = p->nState;
  int nTime = static_cast<int>(std::floor(tmax / dt_out)) + 1;
  Rcpp::NumericMatrix output(nTime, nState + 1);

  std::vector<double> state(x0.begin(), x0.end());

  output(0, 0) = 0.0;
  for (int i = 0; i < nState; i++) output(0, i + 1) = state[i];

  GetRNGstate();

  int out_idx = 1;
  double t = 0.0;
  while (out_idx < nTime) {
    double t_next_out = out_idx * dt_out;
    while (t < t_next_out - 1e-10) {
      double dt = std::min(dt_stoch, t_next_out - t);
      stoch_epi_do_step(p.get(), state.data(), t, dt);
      t += dt;
    }
    output(out_idx, 0) = t;
    for (int i = 0; i < nState; i++) output(out_idx, i + 1) = state[i];
    out_idx++;
  }

  PutRNGstate();
  return output;
}

// ---- Per-step export for R-side simulation framework ----
// Runs sub-steps from t to t + dt_out using dt_stoch steps.
// Returns updated state vector. Brackets with GetRNGstate/PutRNGstate.
// Human state must be set beforehand via stoch_epi_set_human_state() for Imperial.
// [[Rcpp::export]]
Rcpp::NumericVector stoch_epi_step_native(
    SEXP engine_ptr, Rcpp::NumericVector state_r,
    double t, double dt_stoch, double dt_out)
{
  Rcpp::XPtr<StochEpiParams> p(engine_ptr);
  const int nState = p->nState;

  std::vector<double> state(state_r.begin(), state_r.end());

  GetRNGstate();

  double t_now = t;
  double t_end = t + dt_out;

  while (t_now < t_end - 1e-10) {
    double dt = std::min(dt_stoch, t_end - t_now);
    stoch_epi_do_step(p.get(), state.data(), t_now, dt);
    t_now += dt;
  }

  PutRNGstate();

  return Rcpp::NumericVector(state.begin(), state.end());
}

// stoch_engine_ento.cpp
// Native entomology-only stochastic tau-leaping step.

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rmath.h>
using namespace Rcpp;

static inline int ix3s(int stage, int geno, int node, int d1, int d2) {
  return stage + geno * d1 + node * d1 * d2;
}

static inline int fem_row(int gf, int gm, int nG) {
  return gf + gm * nG;
}

static int draw_binom(double n, double p) {
  if (n <= 0.0 || p <= 0.0) return 0;
  if (p > 1.0) p = 1.0;
  return static_cast<int>(R::rbinom(std::floor(n + 0.5), p));
}

static int draw_competing(
    double n, const std::vector<double>& rates, double dt,
    std::vector<int>& out)
{
  std::fill(out.begin(), out.end(), 0);
  if (n <= 0.0) return 0;
  double total_rate = 0.0;
  for (double r : rates) if (r > 0.0) total_rate += r;
  if (total_rate <= 0.0) return 0;

  int n_event = draw_binom(n, 1.0 - std::exp(-total_rate * dt));
  if (n_event <= 0) return 0;

  std::vector<double> probs(rates.size());
  for (size_t i = 0; i < rates.size(); i++)
    probs[i] = (rates[i] > 0.0) ? rates[i] / total_rate : 0.0;
  R::rmultinom(n_event, probs.data(), probs.size(), out.data());
  return n_event;
}

static void alloc_multinom(int n, std::vector<double> probs, std::vector<int>& out) {
  std::fill(out.begin(), out.end(), 0);
  if (n <= 0) return;
  R::rmultinom(n, probs.data(), probs.size(), out.data());
}

static void stoch_ento_one_step(
    std::vector<double>& state,
    int nNodes, int nG, int nE, int nL, int nP,
    const std::vector<int>& egg_ix, const std::vector<int>& larv_ix,
    const std::vector<int>& pup_ix, const std::vector<int>& male_ix,
    const std::vector<int>& unm_ix, const std::vector<int>& fem_ix,
    double rE, double rL, double rP,
    const std::vector<double>& muE, const std::vector<double>& muL,
    const std::vector<double>& muP, const std::vector<double>& muM,
    const std::vector<double>& muF,
    bool log_dd, const std::vector<double>& K, const std::vector<double>& gamma_dd,
    const std::vector<double>& beta_vec, double nu,
    const std::vector<double>& omega_inv, const std::vector<double>& phi,
    const std::vector<double>& xiF, const std::vector<double>& xiM,
    const Rcpp::NumericMatrix& eta, const Rcpp::NumericMatrix& B_mat,
    bool has_move, const Rcpp::NumericMatrix& move_probs,
    const std::vector<double>& move_rates,
    const std::vector<double>& muM_node_base,
    const std::vector<double>& muF_node_base,
    double tol, double dt)
{
  int nPair = nG * nG;
  int nState = (int)state.size();
  for (int i = 0; i < nState; i++)
    if (!R_finite(state[i]) || state[i] < 0.0) state[i] = 0.0;

  std::vector<double> d(state.size(), 0.0);
  std::vector<int> comp2(2), comp3(3), alloc(nG);
  std::vector<double> probs(nG);

  // Births from all mated females.
  for (int k = 0; k < nNodes; k++) {
    for (int go = 0; go < nG; go++) {
      double expected = 0.0;
      for (int gm = 0; gm < nG; gm++) {
        for (int gf = 0; gf < nG; gf++) {
          int pair = fem_row(gf, gm, nG);
          expected += B_mat(pair, go) * state[fem_ix[pair + k * nPair]];
        }
      }
      double lam = beta_vec[k] * expected * dt;
      if (lam > 0.0) d[egg_ix[ix3s(0, go, k, nE, nG)]] += R::rpois(lam);
    }
  }

  // Eggs.
  for (int k = 0; k < nNodes; k++) {
    for (int g = 0; g < nG; g++) {
      for (int e = 0; e < nE; e++) {
        int si = egg_ix[ix3s(e, g, k, nE, nG)];
        draw_competing(state[si], {muE[k], rE}, dt, comp2);
        d[si] -= comp2[0] + comp2[1];
        if (e < nE - 1) {
          d[egg_ix[ix3s(e + 1, g, k, nE, nG)]] += comp2[1];
        } else {
          d[larv_ix[ix3s(0, g, k, nL, nG)]] += comp2[1];
        }
      }
    }
  }

  // Larval density-dependence and progression.
  std::vector<double> dd_rate(nNodes, 0.0);
  for (int k = 0; k < nNodes; k++) {
    double L_tot = 0.0;
    for (int g = 0; g < nG; g++)
      for (int l = 0; l < nL; l++)
        L_tot += state[larv_ix[ix3s(l, g, k, nL, nG)]];
    dd_rate[k] = log_dd ? muL[k] * (1.0 + L_tot / K[k]) : muL[k] + gamma_dd[k] * L_tot;
  }
  for (int k = 0; k < nNodes; k++) {
    for (int g = 0; g < nG; g++) {
      for (int l = 0; l < nL; l++) {
        int si = larv_ix[ix3s(l, g, k, nL, nG)];
        draw_competing(state[si], {dd_rate[k], rL}, dt, comp2);
        d[si] -= comp2[0] + comp2[1];
        if (l < nL - 1) {
          d[larv_ix[ix3s(l + 1, g, k, nL, nG)]] += comp2[1];
        } else {
          d[pup_ix[ix3s(0, g, k, nP, nG)]] += comp2[1];
        }
      }
    }
  }

  // Pupae.
  std::vector<int> male_emerge(nG * nNodes, 0), female_emerge(nG * nNodes, 0);
  for (int k = 0; k < nNodes; k++) {
    double tot_m = 0.0;
    for (int g = 0; g < nG; g++) tot_m += state[male_ix[g + k * nG]];
    for (int g = 0; g < nG; g++) {
      for (int p = 0; p < nP; p++) {
        int si = pup_ix[ix3s(p, g, k, nP, nG)];
        if (p < nP - 1) {
          draw_competing(state[si], {muP[k], rP}, dt, comp2);
          d[si] -= comp2[0] + comp2[1];
          d[pup_ix[ix3s(p + 1, g, k, nP, nG)]] += comp2[1];
        } else {
          double f_rate = rP * phi[g] * xiF[g];
          if (tot_m > tol) {
            double denom = 0.0;
            for (int gm = 0; gm < nG; gm++) denom += eta(g, gm) * state[male_ix[gm + k * nG]];
            if (denom <= tol) f_rate = 0.0;
          }
          double m_rate = rP * (1.0 - phi[g]) * xiM[g];
          draw_competing(state[si], {muP[k], m_rate, f_rate}, dt, comp3);
          d[si] -= comp3[0] + comp3[1] + comp3[2];
          male_emerge[g + k * nG] += comp3[1];
          female_emerge[g + k * nG] += comp3[2];
        }
      }
    }
  }

  // Male deaths and emergence.
  for (int k = 0; k < nNodes; k++) {
    double mu = muM[k] * muM_node_base[k];
    for (int g = 0; g < nG; g++) {
      int si = male_ix[g + k * nG];
      int deaths = draw_binom(state[si], 1.0 - std::exp(-mu * omega_inv[g] * dt));
      d[si] += male_emerge[g + k * nG] - deaths;
    }
  }

  // Females: emergence, unmated death/mating, and mated death.
  for (int k = 0; k < nNodes; k++) {
    double mu = muF[k] * muF_node_base[k];
    for (int gf = 0; gf < nG; gf++) {
      double denom = 0.0;
      for (int gm = 0; gm < nG; gm++) {
        probs[gm] = eta(gf, gm) * state[male_ix[gm + k * nG]];
        denom += probs[gm];
      }
      bool has_mate = denom > tol;
      if (has_mate) {
        for (int gm = 0; gm < nG; gm++) probs[gm] /= denom;
      }

      int emerged = female_emerge[gf + k * nG];
      if (emerged > 0 && has_mate) {
        alloc_multinom(emerged, probs, alloc);
        for (int gm = 0; gm < nG; gm++)
          d[fem_ix[fem_row(gf, gm, nG) + k * nPair]] += alloc[gm];
      } else if (emerged > 0) {
        d[unm_ix[gf + k * nG]] += emerged;
      }

      int ui = unm_ix[gf + k * nG];
      if (has_mate) {
        draw_competing(state[ui], {mu * omega_inv[gf], nu}, dt, comp2);
        d[ui] -= comp2[0] + comp2[1];
        if (comp2[1] > 0) {
          alloc_multinom(comp2[1], probs, alloc);
          for (int gm = 0; gm < nG; gm++)
            d[fem_ix[fem_row(gf, gm, nG) + k * nPair]] += alloc[gm];
        }
      } else {
        int deaths = draw_binom(state[ui], 1.0 - std::exp(-mu * omega_inv[gf] * dt));
        d[ui] -= deaths;
      }

      for (int gm = 0; gm < nG; gm++) {
        int fi = fem_ix[fem_row(gf, gm, nG) + k * nPair];
        int deaths = draw_binom(state[fi], 1.0 - std::exp(-mu * omega_inv[gf] * dt));
        d[fi] -= deaths;
      }
    }
  }

  // Movement: adult males and mated females only.
  if (has_move) {
    for (int origin = 0; origin < nNodes; origin++) {
      double rate = move_rates[origin];
      if (rate <= 0.0) continue;
      std::vector<int> dests;
      std::vector<double> weights;
      double row_sum = 0.0;
      for (int dest = 0; dest < nNodes; dest++) {
        double w = move_probs(origin, dest);
        if (w > 0.0) {
          dests.push_back(dest);
          weights.push_back(w);
          row_sum += w;
        }
      }
      if (dests.empty() || row_sum <= 0.0) continue;
      std::vector<double> move_p(weights.size());
      for (size_t i = 0; i < weights.size(); i++) move_p[i] = weights[i] / row_sum;
      double p_move = 1.0 - std::exp(-rate * row_sum * dt);
      std::vector<int> move_alloc(weights.size());

      for (int g = 0; g < nG; g++) {
        int oi = male_ix[g + origin * nG];
        int movers = draw_binom(state[oi], p_move);
        d[oi] -= movers;
        alloc_multinom(movers, move_p, move_alloc);
        for (size_t i = 0; i < dests.size(); i++)
          d[male_ix[g + dests[i] * nG]] += move_alloc[i];
      }
      for (int pair = 0; pair < nPair; pair++) {
        int oi = fem_ix[pair + origin * nPair];
        int movers = draw_binom(state[oi], p_move);
        d[oi] -= movers;
        alloc_multinom(movers, move_p, move_alloc);
        for (size_t i = 0; i < dests.size(); i++)
          d[fem_ix[pair + dests[i] * nPair]] += move_alloc[i];
      }
    }
  }

  for (int i = 0; i < nState; i++) {
    state[i] += d[i];
    if (!R_finite(state[i]) || state[i] < 0.0) state[i] = 0.0;
    state[i] = std::floor(state[i] + 0.5);
  }
}

// [[Rcpp::export]]
Rcpp::NumericVector stoch_ento_step_native(
    Rcpp::NumericVector state_r, double t, double dt_stoch, double dt_out,
    int nNodes, int nG, int nE, int nL, int nP,
    Rcpp::IntegerVector egg_ix_r, Rcpp::IntegerVector larv_ix_r,
    Rcpp::IntegerVector pup_ix_r, Rcpp::IntegerVector male_ix_r,
    Rcpp::IntegerVector unm_ix_r, Rcpp::IntegerVector fem_ix_r,
    double rE, double rL, double rP,
    Rcpp::NumericVector muE_r, Rcpp::NumericVector muL_r,
    Rcpp::NumericVector muP_r, Rcpp::NumericVector muM_r,
    Rcpp::NumericVector muF_r,
    bool log_dd, Rcpp::NumericVector K_r, Rcpp::NumericVector gamma_dd_r,
    Rcpp::NumericVector beta_vec_r, double nu,
    Rcpp::NumericVector omega_inv_r, Rcpp::NumericVector phi_r,
    Rcpp::NumericVector xiF_r, Rcpp::NumericVector xiM_r,
    Rcpp::NumericMatrix eta_r, Rcpp::NumericMatrix B_mat_r,
    bool has_move, Rcpp::NumericMatrix move_probs_r, Rcpp::NumericVector move_rates_r,
    Rcpp::NumericVector muM_node_base_r, Rcpp::NumericVector muF_node_base_r,
    double tol)
{
  std::vector<double> state(state_r.begin(), state_r.end());
  auto conv = [](Rcpp::IntegerVector v) {
    std::vector<int> out(v.size());
    for (int i = 0; i < v.size(); i++) out[i] = v[i] - 1;
    return out;
  };
  auto to_vec = [](Rcpp::NumericVector v) {
    return std::vector<double>(v.begin(), v.end());
  };

  std::vector<int> egg_ix = conv(egg_ix_r);
  std::vector<int> larv_ix = conv(larv_ix_r);
  std::vector<int> pup_ix = conv(pup_ix_r);
  std::vector<int> male_ix = conv(male_ix_r);
  std::vector<int> unm_ix = conv(unm_ix_r);
  std::vector<int> fem_ix = conv(fem_ix_r);

  std::vector<double> muE = to_vec(muE_r), muL = to_vec(muL_r), muP = to_vec(muP_r);
  std::vector<double> muM = to_vec(muM_r), muF = to_vec(muF_r), K = to_vec(K_r);
  std::vector<double> gamma_dd = to_vec(gamma_dd_r), beta_vec = to_vec(beta_vec_r);
  std::vector<double> omega_inv = to_vec(omega_inv_r), phi = to_vec(phi_r);
  std::vector<double> xiF = to_vec(xiF_r), xiM = to_vec(xiM_r);
  std::vector<double> move_rates = to_vec(move_rates_r);
  std::vector<double> muM_node_base = to_vec(muM_node_base_r);
  std::vector<double> muF_node_base = to_vec(muF_node_base_r);

  double elapsed = 0.0;
  while (elapsed < dt_out - 1e-12) {
    double dt = std::min(dt_stoch, dt_out - elapsed);
    stoch_ento_one_step(
      state, nNodes, nG, nE, nL, nP,
      egg_ix, larv_ix, pup_ix, male_ix, unm_ix, fem_ix,
      rE, rL, rP, muE, muL, muP, muM, muF,
      log_dd, K, gamma_dd, beta_vec, nu, omega_inv, phi, xiF, xiM,
      eta_r, B_mat_r, has_move, move_probs_r, move_rates,
      muM_node_base, muF_node_base, tol, dt);
    elapsed += dt;
  }

  return Rcpp::wrap(state);
}

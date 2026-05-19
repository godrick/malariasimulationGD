/*
 * mosquito_ode.cpp
 *
 *  Created on: 11 Jun 2020
 *      Author: gc1610
 */

#include <Rcpp.h>
#include "aquatic_mosquito_eqs.h"

integration_function_t create_eqs(AquaticMosquitoModel& model) {
  return [&model](const state_t& x, state_t& dxdt, double t) {
    auto kt = model.k_timeseries->at(t, false); 
    auto K = carrying_capacity(
      t,
      model.model_seasonality,
      model.g0,
      model.g,
      model.h,
      kt,
      model.R_bar,
      model.rainfall_floor
    );
    
    auto beta = eggs_laid(model.beta, model.mum, model.f);
    // The compartmental adult solver passes the combined 6-state vector
    // (E, L, P, Sm, Pm, Im) through the aquatic equations. That remains a
    // scalar aquatic system; only the first three entries are aquatic states.
    const bool scalar_aquatic =
      ((x.size() == 3 || x.size() == 6) && model.egg_proportions.size() <= 1);

    if (scalar_aquatic) {
      auto n_larvae = x[get_idx(AquaticState::E)] + x[get_idx(AquaticState::L)];

      if(K == 0){
        // If carrying capacity 0 then all aquatic stages removed
        dxdt[get_idx(AquaticState::E)]  = - x[get_idx(AquaticState::E)];
        dxdt[get_idx(AquaticState::L)] = - x[get_idx(AquaticState::L)];
        dxdt[get_idx(AquaticState::P)] = - x[get_idx(AquaticState::P)];
      } else {
        dxdt[get_idx(AquaticState::E)] = beta * (model.total_M) //new eggs
        - x[get_idx(AquaticState::E)] / model.de //growth to late larval stage
        - x[get_idx(AquaticState::E)] * model.mue * (1 + n_larvae / K); //early larval deaths
        
        dxdt[get_idx(AquaticState::L)] = x[get_idx(AquaticState::E)] / model.de //growth from early larval
        - x[get_idx(AquaticState::L)] / model.dl //growth to pupal
        - x[get_idx(AquaticState::L)] * model.mul * (1 + model.gamma * n_larvae / K); //late larval deaths
        
        dxdt[get_idx(AquaticState::P)] = x[get_idx(AquaticState::L)] / model.dl //growth to pupae
        - x[get_idx(AquaticState::P)] / model.dp //growth to adult
        - x[get_idx(AquaticState::P)] * model.mup; // death of pupae
      }
      return;
    }

    if (x.size() % 3 != 0) {
      Rcpp::stop("Aquatic state length must be a multiple of 3");
    }
    auto G = x.size() / 3;
    std::vector<double> egg_proportions = model.egg_proportions;
    if (egg_proportions.size() != G) {
      if (egg_proportions.size() == 1 && G >= 1) {
        egg_proportions.assign(G, 0.0);
        egg_proportions[0] = 1.0;
      } else {
        Rcpp::stop("egg_proportions length must match genotype-resolved aquatic state dimension");
      }
    }

    double n_larvae = 0.0;
    for (size_t g_idx = 0; g_idx < G; ++g_idx) {
      auto base = 3 * g_idx;
      n_larvae += x[base + get_idx(AquaticState::E)] + x[base + get_idx(AquaticState::L)];
    }

    if(K == 0){
      // If carrying capacity 0 then all aquatic stages removed
      for (size_t g_idx = 0; g_idx < G; ++g_idx) {
        auto base = 3 * g_idx;
        dxdt[base + get_idx(AquaticState::E)] = -x[base + get_idx(AquaticState::E)];
        dxdt[base + get_idx(AquaticState::L)] = -x[base + get_idx(AquaticState::L)];
        dxdt[base + get_idx(AquaticState::P)] = -x[base + get_idx(AquaticState::P)];
      }
    } else {
      auto egg_input_total = beta * model.total_M;
      auto mue_eff_scale = model.mue * (1 + n_larvae / K);
      auto mul_eff_scale = model.mul * (1 + model.gamma * n_larvae / K);
      for (size_t g_idx = 0; g_idx < G; ++g_idx) {
        auto base = 3 * g_idx;
        auto E = x[base + get_idx(AquaticState::E)];
        auto L = x[base + get_idx(AquaticState::L)];
        auto P = x[base + get_idx(AquaticState::P)];
        auto Bg = egg_input_total * egg_proportions[g_idx];

        dxdt[base + get_idx(AquaticState::E)] =
          Bg - E / model.de - E * mue_eff_scale;
        dxdt[base + get_idx(AquaticState::L)] =
          E / model.de - L / model.dl - L * mul_eff_scale;
        dxdt[base + get_idx(AquaticState::P)] =
          L / model.dl - P / model.dp - P * model.mup;
      }
    }
  };
}

AquaticMosquitoModel::AquaticMosquitoModel(
  double beta,
  double de,
  double mue,
  Rcpp::XPtr<Timeseries> k_timeseries,
  double gamma,
  double dl,
  double mul,
  double dp,
  double mup,
  double total_M,
  bool model_seasonality,
  double g0,
  std::vector<double> g,
  std::vector<double> h,
  double R_bar,
  double mum,
  double f,
  double rainfall_floor
):
  beta(beta),
  de(de),
  mue(mue),
  k_timeseries(k_timeseries),
  gamma(gamma),
  dl(dl),
  mul(mul),
  dp(dp),
  mup(mup),
  total_M(total_M),
  model_seasonality(model_seasonality),
  g0(g0),
  g(g),
  h(h),
  R_bar(R_bar),
  mum(mum),
  f(f),
  rainfall_floor(rainfall_floor)
{
  egg_proportions = {1.0};
}



//[[Rcpp::export]]
Rcpp::XPtr<AquaticMosquitoModel> create_aquatic_mosquito_model(
    double beta,
    double de,
    double mue,
    Rcpp::XPtr<Timeseries> k_timeseries,
    double gamma,
    double dl,
    double mul,
    double dp,
    double mup,
    double total_M,
    bool model_seasonality,
    double g0,
    std::vector<double> g,
    std::vector<double> h,
    double R_bar,
    double mum,
    double f,
    double rainfall_floor
) {
  auto model = new AquaticMosquitoModel(
    beta,
    de,
    mue,
    k_timeseries,
    gamma,
    dl,
    mul,
    dp,
    mup,
    total_M,
    model_seasonality,
    g0,
    g,
    h,
    R_bar,
    mum,
    f,
    rainfall_floor
  );
  return Rcpp::XPtr<AquaticMosquitoModel>(model, true);
}

//[[Rcpp::export]]
void aquatic_mosquito_model_update(
    Rcpp::XPtr<AquaticMosquitoModel> model,
    double total_M,
    double f,
    double mum
) {
  model->total_M = total_M;
  model->f = f;
  model->mum = mum;
}

//[[Rcpp::export]]
void aquatic_mosquito_model_set_egg_proportions(
    Rcpp::XPtr<AquaticMosquitoModel> model,
    std::vector<double> egg_proportions
) {
  if (egg_proportions.size() == 0) {
    model->egg_proportions = {1.0};
    return;
  }
  model->egg_proportions = egg_proportions;
}


//[[Rcpp::export]]
Rcpp::XPtr<Solver> create_aquatic_solver(
    Rcpp::XPtr<AquaticMosquitoModel> model,
    std::vector<double> init,
    double r_tol,
    double a_tol,
    size_t max_steps
) {
  return Rcpp::XPtr<Solver>(
    new Solver(init, create_eqs(*model), r_tol, a_tol, max_steps),
    true
  );
}

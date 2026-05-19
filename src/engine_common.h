// engine_common.h
// Shared utilities for native C engines in customMGDrive2
//
// Provides:
// - Seasonal lookup table helpers
// - RK4 integrator for human ODE
// - Common type aliases and index computation utilities

#ifndef ENGINE_COMMON_H
#define ENGINE_COMMON_H

#include <RcppArmadillo.h>
#include <algorithm>  // std::max, std::min
#include <cmath>      // floor

// ---- Seasonal lookup helper -------------------------------------------------
// Clamp time index: t_idx = max(0, min(floor(t), tmax))
inline int seasonal_tidx(double t, int tmax) {
  int idx = static_cast<int>(std::floor(t));
  if (idx < 0) idx = 0;
  if (idx > tmax) idx = tmax;
  return idx;
}

// Lookup scalar parameter from table (length tmax+1)
inline double seasonal_lookup_scalar(const double* tbl, double t, int tmax) {
  return tbl[seasonal_tidx(t, tmax)];
}

// Lookup node-varying parameter from table (nNodes * (tmax+1), column-major)
// Returns pointer to the start of the nNodes-length slice for time t
inline const double* seasonal_lookup_vec(const double* tbl, double t, int tmax, int nNodes) {
  int tidx = seasonal_tidx(t, tmax);
  return tbl + tidx * nNodes;
}

// ---- RK4 integrator ---------------------------------------------------------
// Generic 4th-order Runge-Kutta step for an ODE system.
//
// derivs_fn: function computing dy/dt. Signature:
//   void derivs(int neq, double t, const double* y, double* ydot, void* parms)
//
// neq: number of equations
// t: current time
// dt: step size
// y: input state (length neq) — NOT modified
// y_out: output state (length neq) — may alias y
// work: scratch space (at least 5*neq doubles for rk4_step; 6*neq for rk4_integrate)
// parms: opaque pointer passed to derivs_fn

typedef void (*rk4_derivs_fn)(int neq, double t, const double* y, double* ydot, void* parms);

inline void rk4_step(rk4_derivs_fn derivs, int neq,
                     double t, double dt, const double* y, double* y_out,
                     double* work, void* parms) {
  double* k1 = work;
  double* k2 = work + neq;
  double* k3 = work + 2 * neq;
  double* k4 = work + 3 * neq;

  // We need a temp buffer for intermediate y values; reuse k4 as ytmp after k1 is computed
  // Actually we need ytmp separate from k4. Let's use the space more carefully.
  // work layout: k1[neq], k2[neq], k3[neq], k4[neq] — total 4*neq
  // We need ytmp of size neq for stages 2-4. We can use y_out as ytmp if y_out != y.
  // For safety, require work to be 5*neq and use work+4*neq as ytmp.
  // BUT the doc says 4*neq. Let's just allocate ytmp on top of k4 area after we're done with k3.
  // Actually, the simplest correct approach: require 5*neq work space.

  // Re-documenting: work must be at least 5*neq doubles
  double* ytmp = work + 4 * neq;

  // k1 = f(t, y)
  derivs(neq, t, y, k1, parms);

  // k2 = f(t + dt/2, y + dt/2 * k1)
  double dt2 = dt * 0.5;
  for (int i = 0; i < neq; i++) ytmp[i] = y[i] + dt2 * k1[i];
  derivs(neq, t + dt2, ytmp, k2, parms);

  // k3 = f(t + dt/2, y + dt/2 * k2)
  for (int i = 0; i < neq; i++) ytmp[i] = y[i] + dt2 * k2[i];
  derivs(neq, t + dt2, ytmp, k3, parms);

  // k4 = f(t + dt, y + dt * k3)
  for (int i = 0; i < neq; i++) ytmp[i] = y[i] + dt * k3[i];
  derivs(neq, t + dt, ytmp, k4, parms);

  // y_out = y + dt/6 * (k1 + 2*k2 + 2*k3 + k4)
  double dt6 = dt / 6.0;
  for (int i = 0; i < neq; i++) {
    y_out[i] = y[i] + dt6 * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]);
  }
}

// Multi-step RK4 integration over an interval [t, t+T] with n_steps sub-steps
inline void rk4_integrate(rk4_derivs_fn derivs, int neq,
                          double t0, double T_interval, int n_steps,
                          double* y, double* work, void* parms) {
  double dt = T_interval / n_steps;
  double t = t0;
  // work must be at least 5*neq + neq (for y_tmp swap buffer)
  // Actually rk4_step writes to y_out which can be y itself if we use a temp.
  // Simpler: use a swap buffer at work + 5*neq
  double* y_buf = work + 5 * neq;  // extra neq buffer

  for (int step = 0; step < n_steps; step++) {
    rk4_step(derivs, neq, t, dt, y, y_buf, work, parms);
    // Copy y_buf back to y
    std::copy(y_buf, y_buf + neq, y);
    t += dt;
  }
}

// ---- Index utilities --------------------------------------------------------

// Flatten 3D index (i, j, k) into a linear index for column-major array of dims (d1, d2, d3)
inline int idx3(int i, int j, int k, int d1, int d2) {
  return i + j * d1 + k * d1 * d2;
}

// Flatten 4D index (i, j, k, l) for column-major array of dims (d1, d2, d3, d4)
inline int idx4(int i, int j, int k, int l, int d1, int d2, int d3) {
  return i + j * d1 + k * d1 * d2 + l * d1 * d2 * d3;
}

#endif // ENGINE_COMMON_H

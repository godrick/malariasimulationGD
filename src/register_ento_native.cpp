// register_ento_native.cpp
// Register ento_derivs_native as a C routine for deSolve integration.
// Uses [[Rcpp::init]] so compileAttributes never overwrites this.

#include <Rcpp.h>
#include <R_ext/Rdynload.h>

extern "C" void ento_derivs_native(int*, double*, double*, double*, double*, int*);

// [[Rcpp::init]]
void register_ento_native_routines(DllInfo* dll) {
  static R_NativePrimitiveArgType derivs_types[] = {
    INTSXP, REALSXP, REALSXP, REALSXP, REALSXP, INTSXP
  };

  static R_CMethodDef cMethods[] = {
    {"ento_derivs_native", (DL_FUNC) &ento_derivs_native, 6, derivs_types},
    {NULL, NULL, 0, NULL}
  };

  R_registerRoutines(dll, cMethods, NULL, NULL, NULL);
}

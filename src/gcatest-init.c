#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

SEXP assoc(SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef callMethods[] = {
  {"assoc", (DL_FUNC) &assoc, 4},
  {NULL, NULL, 0}
};

void
R_init_gcatest(DllInfo *info)
{
   R_registerRoutines(info, NULL, callMethods, NULL, NULL);
}

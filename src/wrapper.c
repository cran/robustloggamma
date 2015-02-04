#include <R.h>
#include <Rmath.h>

void F77_SUB(dpgamma)(double *x, double *alph, double *scale, int *lower_tail, int *log_p, double *p) {
  *p = pgamma(*x, *alph, *scale, *lower_tail, *log_p);
}

void F77_SUB(dpnorm5)(double *x, double *mu, double *sigma, int *lower_tail, int *log_p, double *p) {
  *p = pnorm5(*x, *mu, *sigma, *lower_tail, *log_p);
}

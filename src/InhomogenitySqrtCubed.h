#ifndef _INHOMOGENITYSQRTCUBED_H_
#define _INHOMOGENITYSQRTCUBED_H_

void inhomogenity_sqrt_cubed_cv_plus(complex *M_inhom, complex_spline *M_hat,
                                     complex *f, complex *g, double *s,
                                     double s0, double a, double L2, double cut,
                                     int N, int n);

void inhomogenity_sqrt_cubed_cv_minus(complex *M_inhom, complex_spline *M_hat,
                                      complex *f, complex *g, double *s,
                                      double s0, double a, double L2, int N,
                                      int n);

void inhomogenity_sqrt_cubed_below_cut(complex *M_inhom, complex_spline *M_hat,
                                       complex *f, complex *g, double *s,
                                       double s0, double a, double L2, int N,
                                       int n);

void inhomogenity_sqrt_cubed_complex(complex *M_inhom, complex_spline *M_hat,
                                     complex *f, complex *g, complex *s,
                                     char plusminus, double s0, double a,
                                     double L2, int N, int n);

void build_inhomogenity_sqrt_cubed(complex **M_inhom, complex_spline *M_hat,
                                   complex *F, complex *G, double s0, double a,
                                   double L2, double *s_cv_plus,
                                   double *s_cv_minus, double *s_below_cut,
                                   complex **s_cmp, int *N, int n);

#endif

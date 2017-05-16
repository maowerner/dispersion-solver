#ifndef _INHOMOGENITYSQRT_H_
#define _INHOMOGENITYSQRT_H_

void inhomogenity_sqrt_cv_plus(complex *M_inhom, complex_spline *M_hat,
                               complex *f, complex *g, double *s, double s0,
                               double a, double L2, double cut, int N, int n);

void inhomogenity_sqrt_cv_minus(complex *M_inhom, complex_spline *M_hat,
                                complex *f, complex *g, double *s, double s0,
                                double a, double L2, int N, int n);

void inhomogenity_sqrt_below_cut(complex *M_inhom, complex_spline *M_hat,
                                 complex *f, complex *g, double *s, double s0,
                                 double a, double L2, int N, int n);

complex inhomogenity_sqrt_below_cut_single_value(complex_spline *M_hat,
                                                 complex *f, complex *g,
                                                 double s, double s0, double a,
                                                 double L2, int n);

void inhomogenity_sqrt_complex(complex *M_inhom, complex_spline *M_hat,
                               complex *f, complex *g, complex *s,
                               char plusminus, double s0, double a, double L2,
                               int N, int n);

void build_inhomogenity_sqrt(complex **M_inhom, complex_spline *M_hat,
                             complex *F, complex *G, double s0, double a,
                             double L2, double *s_cv_plus, double *s_cv_minus,
                             double *s_below_cut, complex **s_cmp, int *N,
                             int n);

#endif

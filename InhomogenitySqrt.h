#ifndef _INHOMOGENITYSQRT_H_
#define _INHOMOGENITYSQRT_H_

void inhomogenity_sqrt_cv_plus(complex *M_inhom, complex_spline *M_hat, complex *f, complex *g, double *s, complex *R12, complex *Q12, double s0, double L2, double cut, int N, int n);

void inhomogenity_sqrt_cv_minus(complex *M_inhom, complex_spline *M_hat, complex *f, complex *g, double *s, complex *R12, complex *Q12, double s0, double L2, int N, int n);

void inhomogenity_sqrt_below_cut(complex *M_inhom, complex_spline *M_hat, complex *f, complex *g, double *s, complex *Q12, double s0, double L2, int N, int n);

void inhomogenity_sqrt_complex(complex *M_inhom, complex_spline *M_hat, complex *f, complex *g, complex *s, complex *Q12_f, complex *Q12_g, double s0, double L2, int N, int n);

void build_inhomogenity_sqrt(complex **M_inhom, complex_spline *M_hat, complex *F, complex *G, double s0, double L2, double *s_cv_plus, double *s_cv_minus, double *s_below_cut, complex *s_complex, complex *R12_cv_plus, complex *Q12_cv_plus, complex *R12_cv_minus, complex *Q12_cv_minus, complex *Q12_below_cut, complex *Q12_complex_f, complex *Q12_complex_g, int *N, int n);

#endif

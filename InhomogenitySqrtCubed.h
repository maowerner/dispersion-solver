#ifndef _INHOMOGENITYSQRTCUBED_H_
#define _INHOMOGENITYSQRTCUBED_H_

void inhomogenity_sqrt_cubed_cv_plus(complex *M_inhom, complex_spline *M_hat, complex *f, complex *g, double *s, complex *R32, complex *Q12, complex *Q32, double s0, double L2, double cut, int N, int n);

void inhomogenity_sqrt_cubed_cv_minus(complex *M_inhom, complex_spline *M_hat, complex *f, complex *g, double *s, complex *R32, complex *Q12, complex *Q32, double s0, double L2, int N, int n);

void inhomogenity_sqrt_cubed_below_cut(complex *M_inhom, complex_spline *M_hat, complex *f, complex *g, double *s, complex *Q12, complex *Q32, double s0, double L2, int N, int n);

void inhomogenity_sqrt_cubed_complex(complex *M_inhom, complex_spline *M_hat, complex *f, complex *g, complex *s, complex *Q12_f, complex *Q32_f, complex *Q12_g, complex *Q32_g, double s0, double L2, int N, int n);

void build_inhomogenity_sqrt_cubed(complex **M_inhom, complex_spline *M_hat, complex *F, complex *G, double s0, double L2, double *s_cv_plus, double *s_cv_minus, double *s_below_cut, complex *s_complex, complex *R32_cv_plus, complex *Q12_cv_plus, complex *Q32_cv_plus, complex *R32_cv_minus, complex *Q12_cv_minus, complex *Q32_cv_minus, complex *Q12_below_cut, complex *Q32_below_cut, complex *Q12_complex_f, complex *Q32_complex_f, complex *Q12_complex_g, complex *Q32_complex_g, int *N, int n);

#endif

#ifndef _ETAPI_DISC_H_
#define _ETAPI_DISC_H_

void build_inhomogenity_integrand_etapi_disc(complex_spline M0_finite, complex_spline M2_finite, complex_spline *M0_sing, complex_spline *M2_sing, complex *cf0, complex *cg0, complex *cf2, complex *cg2, complex_spline f0etapi, complex_spline f2etapi, complex_spline Metapi, complex_spline *Metapi_tilde, gsl_spline *delta0, gsl_spline *delta2, double *abs_omnes0, double *abs_omnes2, double *s, double s0, double L2, int *N, int Nsin, int n0, int n2);

//etapi_disc_finite function for given s in [s_etapi,L2] above the cut
void etapi_disc_finite_cv_plus(complex *M_etapi, complex_spline M_etapi_hat, double s0, double L2, double *s, int N, int n);

//etapi_disc_finite function for given s in [s_min,s_etapi] below the cut
void etapi_disc_finite_below_cut(complex *M_etapi, complex_spline M_etapi_hat, double s0, double L2, double s_min, double *s, int N, int n);

//etapi_disc_finite function for complex values of s
void etapi_disc_finite_complex(complex *M_etapi, complex_spline M_etapi_hat, double s0, double L2, complex *s, int N, int n);

void build_etapi_inhomogeneity(complex **M0_inhom, complex **M1_inhom, complex **M2_inhom, complex_spline M0_finite, complex_spline M2_finite, complex_spline *M0_sing, complex_spline *M2_sing, complex *F0, complex *F2, complex *G0, complex *G2, double *s_cv_plus, double *s_cv_minus, double *s_below_cut, complex *s_complex, double s_min, double L2, int *N, int n0_sub, int n1_sub, int n2_sub);

#endif

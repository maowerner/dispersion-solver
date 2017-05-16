#ifndef _ANGULARAVERAGES_H_
#define _ANGULARAVERAGES_H_

// lower boundary of the integral s_minus for given s
double integration_s_minus(double s);

// upper boundary of the integral s_plus for given s
double integration_s_plus(double s);

// functions for integration region I, s in [4.0*MPION^2:1/2*(METAP^2-MPION^2)]
// (s' is real +i*eps)

// Mhat in integration region I
void M_avg_1(complex_spline M, complex *M_avg, double *s, double a, double b,
             double c, double d, char plusminus, char kappa_pm, int N, int n);

// functions for integration region II, s in
// [1/2*(METAP^2-MPION^2):(METAP-MPION)^2] (s' is real +i*eps and +i*eps)

// Mhat in integration region II for given s
void M_avg_2(complex_spline M_plus, complex_spline M_minus, complex *M_avg,
             double *s, double a, double b, double c, double d, double s0,
             char plusminus, char kappa_pm, int N, int n);

// functions for integration region III s in [(METAP-MPION)^2:(METAP+MPION)^2]
// (s' is complex)
// parameterising complex path by gamma=r(phi)*exp(i*phi)
// Mhat in integration region III for given s
void M_avg_3(complex_spline M, complex *M_avg, double *s, int N, int n);

void M_avg_III(complex_spline *M, complex *M_avg, double *s, double a, double b,
               double c, double d, char plusminus, char kappa_pm, int N,
               double n);

// functions for integration region IV s in [(METAP+MPION)^2:inf] (s' is real
// and negative)

// Mhat in integration region IV for given s
void M_avg_4(complex_spline M, complex *M_avg, double *s, double a, double b,
             double c, double d, char plusminus, char kappa_pm, int N, int n);

void build_M_avg(complex_spline *M, complex *M_avg, double *s, int *N, double a,
                 double b, double c, double d, double s0, char plusminus,
                 char kappa_pm, int n);

void build_inhomogenity_integrand(complex_spline *M0, complex_spline *M1,
                                  complex_spline *M2, complex *cf0,
                                  complex *cg0, complex *cf1, complex *cg1,
                                  complex *cf2, complex *cg2, complex **M00,
                                  complex **M11, complex **M20, complex *omnes0,
                                  complex *omnes1, complex *omnes2, double *s,
                                  double s0, double L2, int *N, int Nsin,
                                  int n0, int n1, int n2);

void build_inhomogenity_integrand_etap_eta2pi(
    complex_spline *M0, complex_spline *M1, complex *cf0, complex *cg0,
    complex *cf1, complex *cg1, complex *M0_avg, complex *M1_avg_s,
    complex *M1_avg_u, complex *omnes0, complex *omnes1, double *s, double *t,
    double s_th, double s_ps, double a_s, double b_s, double t_th, double t_ps,
    double a_t, double b_t, double L2, int *N_s, int *N_t, int Nsin, int n0,
    int n1);

void build_inhomogenity_integrand_Matrix(
    complex_spline *M0, complex_spline *M1, complex_spline *M2, complex *cf0,
    complex *cg0, complex *cf1, complex *cg1, complex *cf2, complex *cg2,
    complex *M0_tilde, complex *M1_tilde, complex *M2_tilde, double *sin_delta0,
    double *sin_delta1, double *sin_delta2, double *abs_omnes0,
    double *abs_omnes1, double *abs_omnes2, double *s, double s0, double L2,
    int *N, int n0, int n1, int n2);

#endif

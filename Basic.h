#ifndef _BASIC_H_
#define _BASIC_H_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

//Meson Masses in MeV
#define MPION 139.57018
#define MKAON 497.614
#define META 547.862
#define MRHO 775.26
#define METAP 957.78//958.0//547.862//

#define TWO_MPION_SQUARED 4.0//s=(2.0*MPION)^2 in MPION^2
#define LAMBDA_SQUARED 150.0//Integral cutoff

#define DELTA_A 5.e-03 //cutoff for singularity at pseudo-threshold a=METAP_M_MPION_SQUARED

#define TWO_ACOTH_SQRT2 1.762747174039086 //2*acoth(sqrt(2))
#define SQRT2 1.414213562373095 //sqrt(2)

//MKAON^2 in MPION^2
const double MKAON_SQUARED;

//META^2 in MPION^2
const double META_SQUARED;

//MRHO^2 in MPION^2
const double MRHO_SQUARED;

//METAP^2 in MPION^2
const double METAP_SQUARED;

//(METAP-MPION)^2 in MPION^2
const double METAP_M_MPION_SQUARED;

//(METAP+MPION)^2 in MPION^2
const double METAP_P_MPION_SQUARED;

//(META+MPION)^2 in MPION^2
const double s_etapi;

//(META-MPION)^2 in MPION^2
const double s_etapi_pseudo;

//Phase space for etapi-disc
double lambda_etapi(double s);

//string comparison
int str_eq_str(char *str1, char *str2);

//#define SIGMA_MPION( s) (sqrt(1.0-4.0/(s)))
//#define Q2_MPION( s) ((s)/4.0-1.0)

//defines a complex double
typedef  struct {
    double re, im;
} complex;

typedef  struct {
    double r, phi;
} complex_polar;

typedef  struct {
    gsl_spline *re,*im;
} complex_spline;

//sign-function
static inline double sign(double x) {
    if (x>=0.) {
        return 1.;
    }
//    else if (x==0.) {
//        return 0.;
//    }
    else {
        return -1.;
    }
}

//squared sigma_pion function
static inline double SIGMA_SQUARE(double s) {
    return 1.0-4.0/s;
}

//squared CMS momentum for a pion-pion system in MPION^2
static inline double Q_SQUARE(double s) {
    return s/4.0-1.0;
}

//squared CMS momentum for a kaon-kaon system in MPION^2
static inline double Q1_SQUARE(double s) {
    return s/4.0-MKAON_SQUARED;
}

//squared CMS momentum for a eta-eta system in MPION^2
static inline double Q2_SQUARE(double s) {
    return s/4.0-META_SQUARED;
}

//
static inline double W_FUNCTION(double s, double s0) {
    double A,B;
    A = sqrt(s);
    B = sqrt(s0-s);
    return (A-B)/(A+B);
}

complex kappa(double s);

void phase_shifts(gsl_spline *delta0, gsl_spline *delta1, gsl_spline *delta2, double *sin_delta0, double *sin_delta1, double *sin_delta2, double *s, int N);

//angle phi of the complex path for given s
double integration_III_phi_s(double s);

//radius of the complex path for given s
double integration_III_radius_s(double s);

//radius of the complex path for given phi
double integration_III_radius_phi(double phi);

//derivative of the radius with respect to phi for given phi
double integration_III_dradius_phi(double phi);

//path gamma for given phi
complex integration_III_gamma_phi(double phi);

//derivative of the path gamma with respect to phi for given phi
complex integration_III_dgamma_phi(double phi);

double integration_III_circumference();

void complex_spline_alloc(complex_spline *A, int N, int *N_spline);

void complex_spline_free(complex_spline *A, int N);

//void build_s_cv_plus(double *s, double s_init, double s_final, double eps, double eps_b, int *N, int N_sin);
//
//void build_s_cv_minus(double *s, double s_init, double s_final, double eps, int N);
//
//void build_s_below_cut(double *s, double si, double sf, double eps, int N, int N_sin);
//
//void build_s_complex(complex *s, double *phi, double eps, int N);

void cubic_spline_f01_df01(double x0, double x1, double f0, double f1, double df0, double df1, double *c);

complex R_sqrt(double s, double s0, double L2);

complex R_sqrt_cubed(double s, double s0, double L2, complex R12);

complex Q_sqrt(double s, double s0, double L2);

complex Q_sqrt_cubed(double s, double s0, double L2, complex Q12);

//complex Q_sqrt_complex(complex s, double s0, double L2);

complex Q_sqrt_complex_f(complex s, double s0);

complex Q_sqrt_complex_g(complex s, double L2);

complex Q_sqrt_cubed_complex(complex s, double s0, double L2, complex Q12);

complex Q_sqrt_cubed_complex_f(complex s, double s0, complex Q12_f);

complex Q_sqrt_cubed_complex_g(complex s, double L2, complex Q12_g);

#endif

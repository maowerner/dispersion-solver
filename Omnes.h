#ifndef _OMNES_H_
#define _OMNES_H_

complex Omnes_function(gsl_spline *delta, complex s, double s0, double L2, double delta_L2);

//abs(omnes) function for given s in [s0,L2] above the cut for M_hat integration
void abs_omnes_cv_plus(double *abs_omnes, gsl_spline *delta, double s0, double L2, double delta_L2, double *s, int N);

//omnes function for given s in [s0,L2] above the cut
void omnes_cv_plus(complex *omnes, double *abs_omnes, gsl_spline *delta, double s0, double L2, double delta_L2, double *s, int N);

//omnes function for given s in [s0,L2] below the cut
void omnes_cv_minus(complex *omnes, gsl_spline *delta, double s0, double L2, double delta_L2, double *s, int N);

//omnes function for given complex s in complex path for integration region III
void omnes_complex(complex *omnes, gsl_spline *delta, double s0, double L2, double delta_L2, complex *s, int N);

//omnes function for given s in [s_min,s0]
void omnes_below_cut(complex *omnes, gsl_spline *delta, double s0, double L2, double delta_L2, double *s, int N);

//calculating omnes function in all four regions
void build_omnes(complex **omnes, double *abs_omnes, gsl_spline *delta, double s0, double L2, double delta_L2, double *s_cv_plus, double *s_cv_minus, double *s_below_cut, complex *s_complex, int *N);

#endif

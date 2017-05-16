#ifndef _OMNES_H_
#define _OMNES_H_

double dOmnes0(double (*delta)(double), delta_params Delta);

double d2Omnes0(double (*delta)(double), delta_params Delta, double dO0);

complex Omnes_function(double (*delta)(double), delta_params Delta, complex s);

//omnes function for given s in [s0,L2] above the cut
void omnes_cv_plus(complex *omnes, double (*delta)(double), delta_params Delta, double *s, int N);

//omnes function for given s in [s0,L2] below the cut
void omnes_cv_minus(complex *omnes, double (*delta)(double), delta_params Delta, double *s, int N);

//omnes function for given complex s in complex path for integration region III
void omnes_complex(complex *omnes, double (*delta)(double), delta_params Delta, complex *s, int N);

//omnes function for given s in [s_min,s0]
void omnes_below_cut(complex *omnes, double (*delta)(double), delta_params Delta, double *s, int N);

//calculating omnes function in all four regions
void build_omnes(complex **omnes, double (*delta)(double), delta_params Delta, double *s_cv_plus, double *s_cv_minus, double *s_below_cut, complex **s_cmp, int *N);

#endif

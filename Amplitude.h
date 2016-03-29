#ifndef _AMPLITUDE_H_
#define _AMPLITUDE_H_

void subtraction_constant(char *sub_const, int *n0, int *n1, int *n2);

void build_inhomogenity_init(complex **M_inhom, int *N);

void build_amplitude(complex_spline *M, complex **omnes, complex **M_inhom, double *s_cv_plus, double *s_cv_minus, double *s_below_cut, complex *s_complex, double *phi, double s0, double L2, double s_min, int *N, int n);

#endif

#ifndef _AMPLITUDE_H_
#define _AMPLITUDE_H_

void subtraction_constant(char *sub_const, int *n0, int *n1, int *n2);

void build_inhomogenity_init(complex **M_inhom, int *N);

void amplitude_complex(complex_spline M, complex *omnes, complex *M_inhom,
                       complex *s, double *y, char plusminus, int N, int n);

void build_amplitude(complex_spline *M, complex **omnes, complex **M_inhom,
                     double *s_cv_plus, double *s_cv_minus, double *s_below_cut,
                     complex **s_cmp, double **y, double s0, double L2, int *N,
                     int n);

#endif

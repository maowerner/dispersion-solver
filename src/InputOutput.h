#ifndef _INPUTOUTPUT_H_
#define _INPUTOUTPUT_H_

// reading input file
void input(char *filename, double *s_step, double *L2, int *n0_sub, int *n1_sub,
           int *n2_sub);

void phase_input(char *filename, gsl_spline **delta_spline);

void M_tilde_input(char *filename, complex *M0_tilde, complex *M1_tilde,
                   complex *M2_tilde, double *s, int *N_low, int *N_high,
                   double eps_b);

void F_etapi_input(char *filename, double *F0_re, double *F0_im, double *F2_re,
                   double *F2_im, double *s, int *N);

void Partial_wave_eta_3pi_input(char *filename, complex_spline f, double L2,
                                int n);

void Partial_wave_etap_eta2pi_input(char *filename, complex_spline M,
                                    complex_spline *M_tilde, double L2, int n,
                                    int m);

void Output_test(char *filename, complex_spline f0, complex_spline f2,
                 complex_spline M, complex_spline *M_tilde, double ds,
                 double L2);

void Omnes_output(char *filename, complex **omnes0, complex **omnes1,
                  complex **omnes2, double *s_below_cut, double *s_cv_plus,
                  int N_below_cut, int N_cv_plus);

void Inhomogenities_output(char *filename, complex **M0_avg, complex **M1_avg,
                           complex **M2_avg, double *s, int N);

void Dispersion_integral_output(char *filename, complex **M0_inhom,
                                complex **M1_inhom, complex **M2_inhom,
                                double *s_below_cut, double *s_cv_plus,
                                int N_below_cut, int N_cv_plus);

void Amplitudes_output(char *filename, complex_spline *M0, complex_spline *M1,
                       complex_spline *M2, double *s_below_cut,
                       double *s_cv_plus, int N_below_cut, int N_cv_plus);

#endif

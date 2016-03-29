#ifndef _GRID_H_
#define _GRID_H_

void grid_point_size(int *N, int Nsin, double s_step, double s_step_low_high, double s_min, double s_max, double s_th, double eps_b);

void build_grid(double *s_cv_plus, double *s_cv_minus, double *s_below_cut, complex *s_complex, double *phi, double s_min, double s_max, double s_th, double eps_a, double eps_b, int *N, int Nsin);

void build_s_etapi_disc(double *s_etapi_disc, double s_init, double s_final, double eps, int *N, int *N_temp, int Nsin);

#endif

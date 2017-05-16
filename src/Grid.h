#ifndef _GRID_H_
#define _GRID_H_

void grid_point_size(int *N, int Nsin, double s_step, double s_step_low_high,
                     double s_min, double s_max, double s_th, double eps_b);

void Ngrid_etap_eta2pi(int *N, int Nsin, double s_step, double s_step_low_high,
                       double s_min, double s_max, double s_th, double a_s,
                       double b_s, double s_ps, double s_m, double s_12,
                       double t_th, double a_t, double b_t, double t_ps,
                       double t_m, double t_12, double u_m, double u_12,
                       double eps_b);

void grid_etap_eta2pi(double **s, double **ys, double **yt, double **yu,
                      complex **s_cmp, complex **t_cmp, complex **u_cmp, int *N,
                      int Nsin, double s_step, double s_step_low_high,
                      double s_min, double s_max, double s_th, double a_s,
                      double b_s, double s_ps, double s_m, double s_12,
                      double t_th, double a_t, double b_t, double t_ps,
                      double t_m, double t_12, double u_m, double u_12,
                      double eps_b, double eps_a);

void build_grid_etap_eta2pi(double **s, int *N, double *s_below_cut,
                            double *s_cv_plus, double *s_cv_minus,
                            double *t_below_cut, double *t_cv_plus,
                            double *t_cv_minus, bool pi0pi0);

void build_grid(double *s_cv_plus, double *s_cv_minus, double *s_below_cut,
                complex **s_cmp, double **y, double s_min, double s_max,
                double s_th, double eps_a, double eps_b, double a, double b,
                double c, double d, char plusminus, char kappa_pm, int *N,
                int Nsin);

void build_s_etapi_disc(double *s_etapi_disc, double s_init, double s_final,
                        double eps, int *N, int *N_temp, int Nsin);

#endif

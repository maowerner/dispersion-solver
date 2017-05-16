#include "AngularAverages.h"
#include "Basic.h"
#include "PhaseShifts.h"
#include "Singularity.h"

// lower boundary of the integral s_minus for given s
double integration_s_minus(double s) {
  double k;
  k = (1.0 - TWO_MPION_SQUARED / s) * (s - METAP_M_MPION_SQUARED) *
      (s - METAP_P_MPION_SQUARED);
  k = sqrt(k);
  if (s > METAP_P_MPION_SQUARED) {
    k *= -1.0;
  }
  return 0.5 * (METAP_SQUARED + 3.0 - s - k);
}

// upper boundary of the integral s_plus for given s
double integration_s_plus(double s) {
  double k;
  k = (1.0 - TWO_MPION_SQUARED / s) * (s - METAP_M_MPION_SQUARED) *
      (s - METAP_P_MPION_SQUARED);
  k = sqrt(k);
  if (s > METAP_P_MPION_SQUARED) {
    k *= -1.0;
  }
  return 0.5 * (METAP_SQUARED + 3.0 - s + k);
}

// functions for integration region I, s in [4.0*MPION^2:1/2*(METAP^2-MPION^2)]
// (s' is real +i*eps)

// angular average for Mhat in integration region I for given s (<z^n*M> with
// n={0,1})
void M_avg_1(complex_spline M, complex *M_avg, double *s, double a, double b,
             double c, double d, char plusminus, char kappa_pm, int N, int n) {
  int i, N_int;
  double eps, error, s_minus, s_plus, sigma;

  N_int = 1500;
  eps = 1.0e-10;

  gsl_interp_accel *acc = gsl_interp_accel_alloc();

  gsl_integration_workspace *w = gsl_integration_workspace_alloc(N_int);

  gsl_function F;
  gsl_function G;

  double f_re(double z, void *params) {
    double sigma = *(double *)params;
    double f;
    if (n == 0) {
      f = gsl_spline_eval(M.re, z, acc);
    } else if (n == 1) {
      f = (z - sigma) * gsl_spline_eval(M.re, z, acc);
    } else if (n == 2) {
      f = pow(z - sigma, 2.0) * gsl_spline_eval(M.re, z, acc);
    }
    //        if (isnan(f)==1) {
    //            printf("%d %.3e %.3e\n",n,z,sigma);
    //        }
    return f;
  }

  double f_im(double z, void *params) {
    double sigma = *(double *)params;
    double f;
    if (n == 0) {
      f = gsl_spline_eval(M.im, z, acc);
    } else if (n == 1) {
      f = (z - sigma) * gsl_spline_eval(M.im, z, acc);
    } else if (n == 2) {
      f = pow(z - sigma, 2.0) * gsl_spline_eval(M.im, z, acc);
    }
    //        if (isnan(f)==1) {
    //            printf("%d %.3e %.3e\n",n,z,sigma);
    //        }
    return f;
  }

  F.function = &f_re;
  G.function = &f_im;

  for (i = 0; i < N; i++) {
    switch (plusminus) {
    case '+':
      sigma = 0.5 * (0.5 * (a + b + c + d) - s[i] + sqrt(a * b * c * d) / s[i]);
      break;

    case '0':
      sigma = 0.5 * (0.5 * (a + b + c + d) - s[i]);
      break;

    case '-':
      sigma = 0.5 * (0.5 * (a + b + c + d) - s[i] - sqrt(a * b * c * d) / s[i]);
      break;

    default:
      printf("FATAL ERROR: In function M_avg_1, plusminus has to be '+' or '0' "
             "or '-', '%c' is unknown!",
             plusminus);
      exit(1);
      break;
    }

    switch (kappa_pm) {
    case '+':
      s_minus = cmp_re(gamma_path(s[i], a, b, c, d, plusminus, '-'));
      s_plus = cmp_re(gamma_path(s[i], a, b, c, d, plusminus, '+'));
      break;

    case '-':
      s_minus = cmp_re(gamma_path(s[i], a, b, c, d, plusminus, '+'));
      s_plus = cmp_re(gamma_path(s[i], a, b, c, d, plusminus, '-'));
      break;

    default:
      printf("FATAL ERROR: In function M_avg_1, kappa_pm has to be '+' or '-', "
             "'%c' is unknown!",
             kappa_pm);
      break;
    }

    F.params = &sigma;
    G.params = &sigma;

    gsl_integration_qags(&F, s_minus, s_plus, eps, eps, N_int, w, &M_avg[i].re,
                         &error);
    gsl_integration_qags(&G, s_minus, s_plus, eps, eps, N_int, w, &M_avg[i].im,
                         &error);

    if (isnan(M_avg[i].re) == 1) {
      printf("FATAL ERROR: In function M_avg_1: M_avg[i].re=nan, step=%d\n", i);
      printf("sm=%.10e, sp=%.10e, s=%.10e\n", s_minus, s_plus, s[i]);
      exit(1);
    }
    if (isnan(M_avg[i].im) == 1) {
      printf("FATAL ERROR: In function M_avg_1: M_avg[i].im=nan, step=%d\n", i);
      exit(1);
    }
  }

  gsl_integration_workspace_free(w);
  gsl_interp_accel_free(acc);
}

// functions for integration region II, s in
// [1/2*(METAP^2-MPION^2):(METAP-MPION)^2] (s' is real +i*eps and +i*eps)

// Mhat in integration region II for given s
void M_avg_2(complex_spline M_plus, complex_spline M_minus, complex *M_avg,
             double *s, double a, double b, double c, double d, double s0,
             char plusminus, char kappa_pm, int N, int n) {
  int i, N_int;
  double eps, error, s_minus, s_plus, sigma;
  complex temp;

  N_int = 1500;
  eps = 1.0e-10;

  gsl_interp_accel *acc_plus = gsl_interp_accel_alloc();
  gsl_interp_accel *acc_minus = gsl_interp_accel_alloc();

  gsl_integration_workspace *w = gsl_integration_workspace_alloc(N_int);

  gsl_function F;
  gsl_function G;
  gsl_function H;
  gsl_function I;

  double f_re_plus(double x, void *params) {
    double y = *(double *)params;
    double f;
    if (n == 0) {
      f = gsl_spline_eval(M_plus.re, x, acc_plus);
    } else if (n == 1) {
      f = (x - y) * gsl_spline_eval(M_plus.re, x, acc_plus);
    } else if (n == 2) {
      f = pow(x - y, 2.0) * gsl_spline_eval(M_plus.re, x, acc_plus);
    }
    return f;
  }

  double f_im_plus(double x, void *params) {
    double y = *(double *)params;
    double f;
    if (n == 0) {
      f = gsl_spline_eval(M_plus.im, x, acc_plus);
    } else if (n == 1) {
      f = (x - y) * gsl_spline_eval(M_plus.im, x, acc_plus);
    } else if (n == 2) {
      f = pow(x - y, 2.0) * gsl_spline_eval(M_plus.im, x, acc_plus);
    }
    return f;
  }

  double f_re_minus(double x, void *params) {
    double y = *(double *)params;
    double f;
    if (n == 0) {
      f = gsl_spline_eval(M_minus.re, x, acc_minus);
    } else if (n == 1) {
      f = (x - y) * gsl_spline_eval(M_minus.re, x, acc_minus);
    } else if (n == 2) {
      f = pow(x - y, 2.0) * gsl_spline_eval(M_minus.re, x, acc_minus);
    }
    if (isnan(f) == 1) {
      printf("%d %.10e %.10e %.10e\n", n, x, s0,
             gsl_spline_eval(M_minus.re, x, acc_minus));
    }

    return f;
  }

  double f_im_minus(double x, void *params) {
    double y = *(double *)params;
    double f;
    if (n == 0) {
      f = gsl_spline_eval(M_minus.im, x, acc_minus);
    } else if (n == 1) {
      f = (x - y) * gsl_spline_eval(M_minus.im, x, acc_minus);
    } else if (n == 2) {
      f = pow(x - y, 2.0) * gsl_spline_eval(M_minus.im, x, acc_minus);
    }
    if (isnan(f) == 1) {
      printf("%d %.10e %.10e %.10e\n", n, x, s0,
             gsl_spline_eval(M_minus.im, x, acc_minus));
    }
    return f;
  }

  F.function = &f_re_plus;
  G.function = &f_im_plus;
  H.function = &f_re_minus;
  I.function = &f_im_minus;

  for (i = 0; i < N; i++) {
    switch (plusminus) {
    case '+':
      sigma = 0.5 * (0.5 * (a + b + c + d) - s[i] + sqrt(a * b * c * d) / s[i]);
      break;

    case '0':
      sigma = 0.5 * (0.5 * (a + b + c + d) - s[i]);
      break;

    case '-':
      sigma = 0.5 * (0.5 * (a + b + c + d) - s[i] - sqrt(a * b * c * d) / s[i]);
      break;

    default:
      printf("FATAL ERROR: In function M_avg_1, plusminus has to be '+' or '0' "
             "or '-', '%c' is unknown!",
             plusminus);
      exit(1);
      break;
    }

    s_minus = cmp_re(gamma_path(s[i], a, b, c, d, plusminus, '-'));
    s_plus = cmp_re(gamma_path(s[i], a, b, c, d, plusminus, '+'));

    F.params = &sigma;
    gsl_integration_qags(&F, s0 + 1.e-10, s_plus, eps, eps, N_int, w, &temp.re,
                         &error);

    G.params = &sigma;
    gsl_integration_qags(&G, s0 + 1.e-10, s_plus, eps, eps, N_int, w, &temp.im,
                         &error);

    M_avg[i].re = temp.re;
    M_avg[i].im = temp.im;

    H.params = &sigma;
    gsl_integration_qags(&H, s_minus, s0 + 1.e-10, eps, eps, N_int, w, &temp.re,
                         &error);

    I.params = &sigma;
    gsl_integration_qags(&I, s_minus, s0 + 1.e-10, eps, eps, N_int, w, &temp.im,
                         &error);

    M_avg[i].re += temp.re;
    M_avg[i].im += temp.im;

    switch (kappa_pm) {
    case '+':
      break;

    case '-':
      M_avg[i].re *= -1.;
      M_avg[i].im *= -1.;
      break;

    default:
      printf("FATAL ERROR: In function M_avg_1, kappa_pm has to be '+' or '-', "
             "'%c' is unknown!",
             kappa_pm);
      break;
    }

    if (isnan(M_avg[i].re) == 1) {
      printf("FATAL ERROR: In function M_avg_2: M_avg[i].re=nan, step=%d\n", i);
      exit(1);
    }
    if (isnan(M_avg[i].im) == 1) {
      printf("FATAL ERROR: In function M_avg_2: M_avg[i].im=nan, step=%d\n", i);
      // exit(1);
    }
  }

  gsl_integration_workspace_free(w);
  gsl_interp_accel_free(acc_plus);
  gsl_interp_accel_free(acc_minus);
}

// functions for integration region III s in [(METAP-MPION)^2:(METAP+MPION)^2]
// (s' is complex)
// parameterising complex path by gamma=r(phi)*exp(i*phi)

// Mhat in integration region III for given s
void M_avg_3(complex_spline M, complex *M_avg, double *s, int N, int n) {
  int i, N_int;
  double eps, error, phi, sigma;

  N_int = 1500;
  eps = 1.0e-10;

  gsl_interp_accel *acc = gsl_interp_accel_alloc();

  gsl_integration_workspace *w = gsl_integration_workspace_alloc(N_int);

  gsl_function F;
  gsl_function G;

  double f_re(double x, void *params) {
    double y = *(double *)params;
    complex gamma = integration_III_gamma_phi(x);
    complex dgamma = integration_III_dgamma_phi(x);
    double f;
    if (n == 0) {
      f = dgamma.re * gsl_spline_eval(M.re, x, acc) -
          dgamma.im * gsl_spline_eval(M.im, x, acc);
    } else if (n == 1) {
      f = (gamma.re - y) * (dgamma.re * gsl_spline_eval(M.re, x, acc) -
                            dgamma.im * gsl_spline_eval(M.im, x, acc)) -
          gamma.im * (dgamma.re * gsl_spline_eval(M.im, x, acc) +
                      dgamma.im * gsl_spline_eval(M.re, x, acc));
    } else if (n == 2) {
      f = (pow(gamma.re, 2.0) + pow(y, 2.0) - pow(gamma.im, 2.0) -
           2.0 * gamma.re * y) *
              (dgamma.re * gsl_spline_eval(M.re, x, acc) -
               dgamma.im * gsl_spline_eval(M.im, x, acc)) -
          2.0 * gamma.im * (gamma.re - y) *
              (dgamma.re * gsl_spline_eval(M.im, x, acc) +
               dgamma.im * gsl_spline_eval(M.re, x, acc));
    }
    if (isnan(f) == 1) {
      printf("%.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e\n", y, x / M_PI, gamma.re,
             gamma.im, dgamma.re, dgamma.im, gsl_spline_eval(M.re, x, acc),
             gsl_spline_eval(M.im, x, acc));
    }
    return f;
  }

  double f_im(double x, void *params) {
    double y = *(double *)params;
    complex gamma = integration_III_gamma_phi(x);
    complex dgamma = integration_III_dgamma_phi(x);
    double f;
    if (n == 0) {
      f = dgamma.re * gsl_spline_eval(M.im, x, acc) +
          dgamma.im * gsl_spline_eval(M.re, x, acc);
    } else if (n == 1) {
      f = (gamma.re - y) * (dgamma.re * gsl_spline_eval(M.im, x, acc) +
                            dgamma.im * gsl_spline_eval(M.re, x, acc)) +
          gamma.im * (dgamma.re * gsl_spline_eval(M.re, x, acc) -
                      dgamma.im * gsl_spline_eval(M.im, x, acc));
    } else if (n == 2) {
      f = (pow(gamma.re, 2.0) + pow(y, 2.0) - pow(gamma.im, 2.0) -
           2.0 * gamma.re * y) *
              (dgamma.re * gsl_spline_eval(M.im, x, acc) +
               dgamma.im * gsl_spline_eval(M.re, x, acc)) +
          2.0 * gamma.im * (gamma.re - y) *
              (dgamma.re * gsl_spline_eval(M.re, x, acc) -
               dgamma.im * gsl_spline_eval(M.im, x, acc));
    }
    return f;
  }

  F.function = &f_re;
  G.function = &f_im;

  for (i = 0; i < N; i++) {
    phi = integration_III_phi_s(s[i]);

    sigma = 0.5 * (METAP_SQUARED + 3.0 - s[i]);

    F.params = &sigma;
    gsl_integration_qags(&F, 2.0 * M_PI - phi, phi, eps, eps, N_int, w,
                         &M_avg[i].re, &error);

    G.params = &sigma;
    gsl_integration_qags(&G, 2.0 * M_PI - phi, phi, eps, eps, N_int, w,
                         &M_avg[i].im, &error);

    if (isnan(M_avg[i].re) == 1) {
      printf("FATAL ERROR: In function M_avg_3: M_avg[i].re=nan, step=%d\n", i);
      exit(1);
    }
    if (isnan(M_avg[i].im) == 1) {
      printf("FATAL ERROR: In function M_avg_3: M_avg[i].im=nan, step=%d\n", i);
      exit(1);
    }
  }

  gsl_integration_workspace_free(w);
  gsl_interp_accel_free(acc);
}

complex gamma_param(double s) {
  complex gamma;

  gamma.re = 0.5 * (METAP_SQUARED + 3.0 - s);
  gamma.im = 0.5 * sqrt(1. - 4. / s) *
             sqrt((METAP_M_MPION_SQUARED - s) * (s - METAP_P_MPION_SQUARED));

  return gamma;
}

complex dgamma_param(double s) {
  complex dgamma;

  dgamma.re = -0.5;
  dgamma.im =
      (4. * METAP_M_MPION_SQUARED * METAP_P_MPION_SQUARED -
       (4. + METAP_M_MPION_SQUARED + METAP_P_MPION_SQUARED) * s * s +
       2. * s * s * s) /
      (4. * sqrt((METAP_M_MPION_SQUARED - s) * (s - METAP_P_MPION_SQUARED)) *
       sqrt(s - 4.) * pow(s, 1.5));

  return dgamma;
}

void M_avg_III(complex_spline *M, complex *M_avg, double *s, double a, double b,
               double c, double d, char plusminus, char kappa_pm, int N,
               double n) {
  int i, N_int;
  double eps, error, sigma, temp, s_ab_avg, yb, yab_avg, yba_avg, yi, y, r, phi,
      c1, s1, c2, s2, f;

  complex gamma, dgamma, M_sub, M_add;

  s_ab_avg = 0.5 * (a + b);

  yb = pow(b, 1. / 4.);
  yab_avg = sqrt(0.5 * (sqrt(s_ab_avg) + sqrt(a))) -
            sqrt(0.5 * (sqrt(s_ab_avg) - sqrt(a)));
  yba_avg = sqrt(0.5 * (sqrt(b) + sqrt(s_ab_avg))) -
            sqrt(0.5 * (sqrt(b) - sqrt(s_ab_avg)));

  N_int = 1500;
  eps = 1.0e-10;

  gsl_interp_accel *acc_a = gsl_interp_accel_alloc();
  gsl_interp_accel *acc_b = gsl_interp_accel_alloc();

  gsl_integration_workspace *w = gsl_integration_workspace_alloc(N_int);

  gsl_function Fa_re;
  gsl_function Fa_im;

  gsl_function Fb_re;
  gsl_function Fb_im;

  double fa_re(double x, void *params) {
    y = *(double *)params;
    gamma = gamma_III(x, a, b, c, d, 'a', plusminus, kappa_pm);
    dgamma = dgamma_III(x, a, b, c, d, 'a', plusminus, kappa_pm);

    r = sqrt(pow(gamma.re - y, 2.) + gamma.im * gamma.im);

    phi = atan2(-gamma.im, gamma.re - y);
    c1 = cos(n * phi);
    s1 = sin(n * phi);

    phi = atan2(gamma.im, gamma.re - y);
    c2 = cos(n * phi);
    s2 = sin(n * phi);

    M_sub.re = c1 * gsl_spline_eval(M[0].re, x, acc_a) -
               c2 * gsl_spline_eval(M[2].re, x, acc_a);
    M_sub.im = s1 * gsl_spline_eval(M[0].im, x, acc_a) -
               s2 * gsl_spline_eval(M[2].im, x, acc_a);

    M_add.re = s1 * gsl_spline_eval(M[0].re, x, acc_a) +
               s2 * gsl_spline_eval(M[2].re, x, acc_a);
    M_add.im = c1 * gsl_spline_eval(M[0].im, x, acc_a) +
               c2 * gsl_spline_eval(M[2].im, x, acc_a);

    f = pow(r, n) * (dgamma.re * M_sub.re + dgamma.im * M_add.im -
                     dgamma.re * M_sub.im + dgamma.im * M_add.re);

    if (isnan(f) == 1) {
      printf("fa_re: %.6e %.6e %.6e %.6e %.6e\n", x, dgamma.re, dgamma.im,
             M_sub.re, M_sub.im);
    }

    return f;
  }

  double fb_re(double x, void *params) {
    y = *(double *)params;
    gamma = gamma_III(x, a, b, c, d, 'b', plusminus, kappa_pm);
    dgamma = dgamma_III(x, a, b, c, d, 'b', plusminus, kappa_pm);

    r = sqrt(pow(gamma.re - y, 2.) + gamma.im * gamma.im);

    phi = atan2(-gamma.im, gamma.re - y);
    c1 = cos(n * phi);
    s1 = sin(n * phi);

    phi = atan2(gamma.im, gamma.re - y);
    c2 = cos(n * phi);
    s2 = sin(n * phi);

    M_sub.re = c1 * gsl_spline_eval(M[1].re, x, acc_b) -
               c2 * gsl_spline_eval(M[3].re, x, acc_b);
    M_sub.im = s1 * gsl_spline_eval(M[1].im, x, acc_b) -
               s2 * gsl_spline_eval(M[3].im, x, acc_b);

    M_add.re = s1 * gsl_spline_eval(M[1].re, x, acc_b) +
               s2 * gsl_spline_eval(M[3].re, x, acc_b);
    M_add.im = c1 * gsl_spline_eval(M[1].im, x, acc_b) +
               c2 * gsl_spline_eval(M[3].im, x, acc_b);

    f = pow(r, n) * (dgamma.re * M_sub.re + dgamma.im * M_add.im -
                     dgamma.re * M_sub.im + dgamma.im * M_add.re);

    if (isnan(f) == 1) {
      printf("fb_re: %.6e %.6e %.6e %.6e %.6e\n", x, dgamma.re, dgamma.im,
             M_sub.re, M_sub.im);
    }

    return f;
  }

  double fa_im(double x, void *params) {
    y = *(double *)params;
    gamma = gamma_III(x, a, b, c, d, 'a', plusminus, kappa_pm);
    dgamma = dgamma_III(x, a, b, c, d, 'a', plusminus, kappa_pm);

    r = sqrt(pow(gamma.re - y, 2.) + gamma.im * gamma.im);

    phi = atan2(-gamma.im, gamma.re - y);
    c1 = cos(n * phi);
    s1 = sin(n * phi);

    phi = atan2(gamma.im, gamma.re - y);
    c2 = cos(n * phi);
    s2 = sin(n * phi);

    M_sub.re = s1 * gsl_spline_eval(M[0].re, x, acc_a) -
               s2 * gsl_spline_eval(M[2].re, x, acc_a);
    M_sub.im = c1 * gsl_spline_eval(M[0].im, x, acc_a) -
               c2 * gsl_spline_eval(M[2].im, x, acc_a);

    M_add.re = c1 * gsl_spline_eval(M[0].re, x, acc_a) +
               c2 * gsl_spline_eval(M[2].re, x, acc_a);
    M_add.im = s1 * gsl_spline_eval(M[0].im, x, acc_a) +
               s2 * gsl_spline_eval(M[2].im, x, acc_a);

    f = pow(r, n) * (dgamma.re * M_sub.im - dgamma.im * M_add.re +
                     dgamma.re * M_sub.re + dgamma.im * M_add.im);

    if (isnan(f) == 1) {
      printf("fa_im: %.6e %.6e %.6e %.6e %.6e\n", x, dgamma.re, dgamma.im,
             M_sub.re, M_sub.im);
    }

    return f;
  }

  double fb_im(double x, void *params) {
    y = *(double *)params;
    gamma = gamma_III(x, a, b, c, d, 'b', plusminus, kappa_pm);
    dgamma = dgamma_III(x, a, b, c, d, 'b', plusminus, kappa_pm);

    r = sqrt(pow(gamma.re - y, 2.) + gamma.im * gamma.im);

    phi = atan2(-gamma.im, gamma.re - y);
    c1 = cos(n * phi);
    s1 = sin(n * phi);

    phi = atan2(gamma.im, gamma.re - y);
    c2 = cos(n * phi);
    s2 = sin(n * phi);

    M_sub.re = s1 * gsl_spline_eval(M[1].re, x, acc_b) -
               s2 * gsl_spline_eval(M[3].re, x, acc_b);
    M_sub.im = c1 * gsl_spline_eval(M[1].im, x, acc_b) -
               c2 * gsl_spline_eval(M[3].im, x, acc_b);

    M_add.re = c1 * gsl_spline_eval(M[1].re, x, acc_b) +
               c2 * gsl_spline_eval(M[3].re, x, acc_b);
    M_add.im = s1 * gsl_spline_eval(M[1].im, x, acc_b) +
               s2 * gsl_spline_eval(M[3].im, x, acc_b);

    f = pow(r, n) * (dgamma.re * M_sub.im - dgamma.im * M_add.re +
                     dgamma.re * M_sub.re + dgamma.im * M_add.im);

    if (isnan(f) == 1) {
      printf("fb_im: %.6e %.6e %.6e %.6e %.6e\n", x, dgamma.re, dgamma.im,
             M_sub.re, M_sub.im);
    }

    return f;
  }

  Fa_re.function = &fa_re;
  Fa_im.function = &fa_im;
  Fb_re.function = &fb_re;
  Fb_im.function = &fb_im;

  for (i = 0; i < N; i++) {
    switch (plusminus) {
    case '+':
      sigma = 0.5 * (0.5 * (a + b + c + d) - s[i] + sqrt(a * b * c * d) / s[i]);
      break;

    case '0':
      sigma = 0.5 * (0.5 * (a + b + c + d) - s[i]);
      break;

    case '-':
      sigma = 0.5 * (0.5 * (a + b + c + d) - s[i] - sqrt(a * b * c * d) / s[i]);
      break;

    default:
      printf("FATAL ERROR: In function M_avg_III, plusminus has to be '+' or "
             "'0' or '-', '%c' is unknown!",
             plusminus);
      exit(1);
      break;
    }

    Fa_re.params = &sigma;
    Fa_im.params = &sigma;
    Fb_re.params = &sigma;
    Fb_im.params = &sigma;

    if (s[i] > s_ab_avg) {
      yi = s_to_y(s[i], a, b);

      gsl_integration_qags(&Fb_re, yi, yb, eps, eps, N_int, w, &M_avg[i].re,
                           &error);
      gsl_integration_qags(&Fb_im, yi, yb, eps, eps, N_int, w, &M_avg[i].im,
                           &error);
    }

    else {
      yi = s_to_y(s[i], a, b);

      gsl_integration_qags(&Fb_re, yba_avg, yb, eps, eps, N_int, w,
                           &M_avg[i].re, &error);
      gsl_integration_qags(&Fb_im, yba_avg, yb, eps, eps, N_int, w,
                           &M_avg[i].im, &error);

      gsl_integration_qags(&Fa_re, yi, yab_avg, eps, eps, N_int, w, &temp,
                           &error);
      M_avg[i].re += temp;
      gsl_integration_qags(&Fa_im, yi, yab_avg, eps, eps, N_int, w, &temp,
                           &error);
      M_avg[i].im += temp;
    }

    if (isnan(M_avg[i].re) == 1) {
      printf("FATAL ERROR: In function M_avg_III: M_avg[i].re=nan, step=%d\n",
             i);
      exit(1);
    }
    if (isnan(M_avg[i].im) == 1) {
      printf("FATAL ERROR: In function M_avg_III: M_avg[i].im=nan, step=%d\n",
             i);
      exit(1);
    }
  }

  gsl_integration_workspace_free(w);
  gsl_interp_accel_free(acc_a);
  gsl_interp_accel_free(acc_b);
}

// functions for integration region IV s in [(METAP+MPION)^2:inf] (s' is real
// and negative)

// Mhat in integration region IV for given s
void M_avg_4(complex_spline M, complex *M_avg, double *s, double a, double b,
             double c, double d, char plusminus, char kappa_pm, int N, int n) {
  int i, N_int;
  double eps, error, s_minus, s_plus, sigma, result;

  N_int = 1500;
  eps = 1.0e-10;

  gsl_interp_accel *acc = gsl_interp_accel_alloc();

  gsl_integration_workspace *w = gsl_integration_workspace_alloc(N_int);

  gsl_function F;
  gsl_function G;

  double f_re(double z, void *params) {
    double sigma = *(double *)params;
    double f;
    if (n == 0) {
      f = gsl_spline_eval(M.re, z, acc);
    } else if (n == 1) {
      f = (z - sigma) * gsl_spline_eval(M.re, z, acc);
    } else if (n == 2) {
      f = pow(z - sigma, 2.0) * gsl_spline_eval(M.re, z, acc);
    }

    if (isnan(f) == 1) {
      printf("f_re: %.6e %.6e\n", z, gsl_spline_eval(M.re, z, acc));
    }

    return f;
  }

  double f_im(double z, void *params) {
    double sigma = *(double *)params;
    double f;
    if (n == 0) {
      f = gsl_spline_eval(M.im, z, acc);
    } else if (n == 1) {
      f = (z - sigma) * gsl_spline_eval(M.im, z, acc);
    } else if (n == 2) {
      f = pow(z - sigma, 2.0) * gsl_spline_eval(M.im, z, acc);
    }

    if (isnan(f) == 1) {
      printf("f_im: %.6e %.6e\n", z, gsl_spline_eval(M.im, z, acc));
    }

    return f;
  }

  F.function = &f_re;
  G.function = &f_im;

  for (i = 0; i < N; i++) {
    switch (plusminus) {
    case '+':
      sigma = 0.5 * (0.5 * (a + b + c + d) - s[i] + sqrt(a * b * c * d) / s[i]);
      break;

    case '0':
      sigma = 0.5 * (0.5 * (a + b + c + d) - s[i]);
      break;

    case '-':
      sigma = 0.5 * (0.5 * (a + b + c + d) - s[i] - sqrt(a * b * c * d) / s[i]);
      break;

    default:
      printf("FATAL ERROR: In function M_avg_1, plusminus has to be '+' or '0' "
             "or '-', '%c' is unknown!",
             plusminus);
      exit(1);
      break;
    }

    switch (kappa_pm) {
    case '+':
      s_minus = cmp_re(gamma_path(s[i], a, b, c, d, plusminus, '-'));
      s_plus = cmp_re(gamma_path(s[i], a, b, c, d, plusminus, '+'));
      break;

    case '-':
      s_minus = cmp_re(gamma_path(s[i], a, b, c, d, plusminus, '+'));
      s_plus = cmp_re(gamma_path(s[i], a, b, c, d, plusminus, '-'));
      break;

    default:
      printf("FATAL ERROR: In function M_avg_1, kappa_pm has to be '+' or '-', "
             "'%c' is unknown!",
             kappa_pm);
      break;
    }

    F.params = &sigma;
    G.params = &sigma;

    gsl_integration_qags(&F, s_minus, s_plus, eps, eps, N_int, w, &M_avg[i].re,
                         &error);
    gsl_integration_qags(&G, s_minus, s_plus, eps, eps, N_int, w, &M_avg[i].im,
                         &error);

    if (isnan(M_avg[i].re) == 1) {
      printf("FATAL ERROR: In function M_avg_4: M_avg[i].re=nan, step=%d\n", i);
      exit(1);
    }
    if (isnan(M_avg[i].im) == 1) {
      printf("FATAL ERROR: In function M_avg_4: M_avg[i].im=nan, step=%d\n", i);
      exit(1);
    }
  }

  gsl_integration_workspace_free(w);
  gsl_interp_accel_free(acc);
}

void build_M_avg(complex_spline *M, complex *M_avg, double *s, int *N, double a,
                 double b, double c, double d, double s0, char plusminus,
                 char kappa_pm, int n) {
  M_avg_1(M[0], M_avg, s, a, b, c, d, plusminus, kappa_pm, N[0], n);
  M_avg_2(M[0], M[1], M_avg + N[0], s + N[0], a, b, c, d, s0, plusminus,
          kappa_pm, N[1], n);
  M_avg_III(M + 3, M_avg + N[0] + N[1], s + N[0] + N[1], a, b, c, d, plusminus,
            kappa_pm, N[2], n);
  M_avg_4(M[2], M_avg + N[0] + N[1] + N[2], s + N[0] + N[1] + N[2], a, b, c, d,
          plusminus, kappa_pm, N[3], n);
}

//
void build_inhomogenity_integrand(complex_spline *M0, complex_spline *M1,
                                  complex_spline *M2, complex *cf0,
                                  complex *cg0, complex *cf1, complex *cg1,
                                  complex *cf2, complex *cg2, complex **M00,
                                  complex **M11, complex **M20, complex *omnes0,
                                  complex *omnes1, complex *omnes2, double *s,
                                  double s0, double L2, int *N, int Nsin,
                                  int n0, int n1, int n2) {
  int i, j, n, m, N1, N2;
  double si, sf, s_am, s_ap, stu, k, c0, c1, c2, *M0_re, *M0_im, *M1_re, *M1_im,
      *M2_re, *M2_im, s_f[5], s_g[5];
  complex temp0, temp1, temp2, f0[4], g0[4], f1[5], g1[5], f2[4], g2[4];

  M0_re = (double *)malloc((N[0] + N[1] + N[2] + N[3]) * sizeof(double));
  M0_im = (double *)malloc((N[0] + N[1] + N[2] + N[3]) * sizeof(double));

  M1_re = (double *)malloc((N[0] + N[1] + N[2] + N[3]) * sizeof(double));
  M1_im = (double *)malloc((N[0] + N[1] + N[2] + N[3]) * sizeof(double));

  M2_re = (double *)malloc((N[0] + N[1] + N[2] + N[3]) * sizeof(double));
  M2_im = (double *)malloc((N[0] + N[1] + N[2] + N[3]) * sizeof(double));

  n = 0;
  for (j = 0; j < 4; j++) {
    if (j == 0) {
      n = 0;
    } else {
      n += N[j - 1];
    }

    m = 0;
    for (i = 0; i < N[j]; i++) {

      stu = s[i + n] - (METAP_SQUARED / 3.0 + 1.0);

      temp0.re = 2.0 / 3.0 * M00[0][i + n].re + 2.0 * stu * M11[0][i + n].re +
                 4.0 / 3.0 * M11[1][i + n].re + 20.0 / 9.0 * M20[0][i + n].re;
      temp0.im = 2.0 / 3.0 * M00[0][i + n].im + 2.0 * stu * M11[0][i + n].im +
                 4.0 / 3.0 * M11[1][i + n].im + 20.0 / 9.0 * M20[0][i + n].im;

      temp1.re = 6.0 * M00[1][i + n].re + 9.0 * stu * M11[1][i + n].re +
                 6.0 * M11[2][i + n].re - 10.0 * M20[1][i + n].re;
      temp1.im = 6.0 * M00[1][i + n].im + 9.0 * stu * M11[1][i + n].im +
                 6.0 * M11[2][i + n].im - 10.0 * M20[1][i + n].im;

      temp2.re = M00[0][i + n].re - 3.0 / 2.0 * stu * M11[0][i + n].re -
                 M11[1][i + n].re + 1.0 / 3.0 * M20[0][i + n].re;
      temp2.im = M00[0][i + n].im - 3.0 / 2.0 * stu * M11[0][i + n].im -
                 M11[1][i + n].im + 1.0 / 3.0 * M20[0][i + n].im;

      k = sqrt(1.0 - 4.0 / s[i + n]) *
          sqrt(fabs(s[i + n] - METAP_P_MPION_SQUARED));

      c0 = sin(delta0(s[i + n])) / (pow(s[i + n], (double)n0) * k *
                                    sqrt(omnes0[i + n].re * omnes0[i + n].re +
                                         omnes0[i + n].im * omnes0[i + n].im));
      c1 = sin(delta1(s[i + n])) / (pow(s[i + n], (double)n1) * k * k * k *
                                    sqrt(omnes1[i + n].re * omnes1[i + n].re +
                                         omnes1[i + n].im * omnes1[i + n].im));
      c2 = sin(delta2(s[i + n])) / (pow(s[i + n], (double)n2) * k *
                                    sqrt(omnes2[i + n].re * omnes2[i + n].re +
                                         omnes2[i + n].im * omnes2[i + n].im));

      if (j == 0 || j == 1) {
        M0_re[i + n] = c0 * temp0.re;
        M0_im[i + n] = c0 * temp0.im;
        M1_re[i + n] = c1 * temp1.re;
        M1_im[i + n] = c1 * temp1.im;
        M2_re[i + n] = c2 * temp2.re;
        M2_im[i + n] = c2 * temp2.im;
      }

      else if (j == 2) {
        // printf("%.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e
        // %.10e %.10e %.10e %.10e
        // %.10e\n",s[i+n],M00[0][i+n].re,M00[0][i+n].im,M00[1][i+n].re,M00[1][i+n].im,M11[0][i+n].re,M11[0][i+n].im,M11[1][i+n].re,M11[1][i+n].im,M11[2][i+n].re,M11[2][i+n].im,M20[0][i+n].re,M20[0][i+n].im,M20[1][i+n].re,M20[1][i+n].im);
        M0_re[i + n] = c0 * temp0.im;
        M0_im[i + n] = -c0 * temp0.re;
        M1_re[i + n] = -c1 * temp1.im;
        M1_im[i + n] = c1 * temp1.re;
        M2_re[i + n] = c2 * temp2.im;
        M2_im[i + n] = -c2 * temp2.re;
      }

      else if (j == 3) {
        M0_re[i + n] = -c0 * temp0.re;
        M0_im[i + n] = -c0 * temp0.im;
        M1_re[i + n] = -c1 * temp1.re;
        M1_im[i + n] = -c1 * temp1.im;
        M2_re[i + n] = -c2 * temp2.re;
        M2_im[i + n] = -c2 * temp2.im;
      }

      // printf("%.10e %.10e %.10e %.10e
      // %.10e\n",s[i+n],M00[0][i+n].re,M00[0][i+n].im,M00[1][i+n].re,M00[1][i+n].im);

      // printf("%+.10e %+.10e %+.10e %+.10e %+.10e %+.10e
      // %+.10e\n",s[i+n],M0_re[i+n],M0_im[i+n],M1_re[i+n],M1_im[i+n],M2_re[i+n],M2_im[i+n]);
      // printf("%+.10e %+.10e %+.10e %+.10e %+.10e %+.10e
      // %+.10e\n",s[i+n],temp0.re,temp0.im,temp1.re,temp1.im,temp2.re,temp2.im);
      // printf("%+.10e %+.10e %+.10e
      // %+.10e\n",s[i+n],sin_delta0[i+n],sin_delta1[i+n],sin_delta2[i+n]);
      // printf("%+.10e %+.10e %+.10e
      // %+.10e\n",s[i+n],abs_omnes0[i+n],abs_omnes1[i+n],abs_omnes2[i+n]);
      if (isnan(temp0.re) == 1) {
        printf("FATAL ERROR: M0_hat.re=nan, step=%d\n", i);
        exit(1);
      }
      if (isnan(temp0.im) == 1) {
        printf("FATAL ERROR: M0_hat.im=nan, step=%d\n", i);
        exit(1);
      }
      if (isnan(temp1.re) == 1) {
        printf("FATAL ERROR: M1_hat.re=nan, step=%d\n", i);
        exit(1);
      }
      if (isnan(temp1.im) == 1) {
        printf("FATAL ERROR: M1_hat.im=nan, step=%d\n", i);
        exit(1);
      }
      if (isnan(temp2.re) == 1) {
        printf("FATAL ERROR: M2_hat.re=nan, step=%d\n", i);
        exit(1);
      }
      if (isnan(temp2.im) == 1) {
        printf("FATAL ERROR: M2_hat.im=nan, step=%d\n", i);
        exit(1);
      }
    }
  }

  m = 0;
  N1 = N[0] + N[1] - 1;
  N2 = N[0] + N[1];
  for (i = 0; i < 4; i++) {
    f0[i].re = M0_re[N1 - m];
    f0[i].im = M0_im[N1 - m];

    f1[i].re = M1_re[N1 - m];
    f1[i].im = M1_im[N1 - m];

    f2[i].re = M2_re[N1 - m];
    f2[i].im = M2_im[N1 - m];

    g0[i].re = M0_re[N2 + i];
    g0[i].im = M0_im[N2 + i];

    g1[i].re = M1_re[N2 + i];
    g1[i].im = M1_im[N2 + i];

    g2[i].re = M2_re[N2 + i];
    g2[i].im = M2_im[N2 + i];

    s_f[i] = sqrt(fabs(s[N1 - m] - METAP_M_MPION_SQUARED));
    s_g[i] = sqrt(fabs(s[N2 + i] - METAP_M_MPION_SQUARED));
    m++;
    if (i == 0) {
      N1 -= Nsin;
      N2 += Nsin;
    }
  }
  f1[i].re = M1_re[N1 - m];
  f1[i].im = M1_im[N1 - m];

  g1[i].re = M1_re[N2 + i];
  g1[i].im = M1_im[N2 + i];

  s_f[i] = sqrt(fabs(s[N1 - m] - METAP_M_MPION_SQUARED));
  s_g[i] = sqrt(fabs(s[N2 + i] - METAP_M_MPION_SQUARED));

  singularity_sqrt(f0, s_f, cf0);
  singularity_sqrt(g0, s_g, cg0);

  // singularity_sqrt(f1,s_f,cf1);
  // singularity_sqrt(g1,s_g,cg1);
  singularity_sqrt_cubed(f1, s_f, cf1);
  singularity_sqrt_cubed(g1, s_g, cg1);

  singularity_sqrt(f2, s_f, cf2);
  singularity_sqrt(g2, s_g, cg2);

  // printf("fre(x)=%+.5e%+.5e*(a-x)%+.5e*(a-x)**1.5%+.5e*(a-x)**2.0%+.5e*(a-x)**2.5\n",cf1[0].re,cf1[1].re,cf1[2].re,cf1[3].re,cf1[4].re);
  // printf("fim(x)=%+.5e%+.5e*(a-x)%+.5e*(a-x)**1.5%+.5e*(a-x)**2.0%+.5e*(a-x)**2.5\n",cf1[0].im,cf1[1].im,cf1[2].im,cf1[3].im,cf1[4].im);
  // printf("gre(x)=%+.5e%+.5e*(x-a)%+.5e*(x-a)**1.5%+.5e*(x-a)**2.0%+.5e*(x-a)**2.5\n",cg1[0].re,cg1[1].re,cg1[2].re,cg1[3].re,cg1[4].re);
  // printf("gim(x)=%+.5e%+.5e*(x-a)%+.5e*(x-a)**1.5%+.5e*(x-a)**2.0%+.5e*(x-a)**2.5\n",cg1[0].im,cg1[1].im,cg1[2].im,cg1[3].im,cg1[4].im);

  // printf("f(x)=%+.5e%+.5e*sqrt(a-x)%+.5e*(a-x)%+.5e*(a-x)**1.5\n",cf0[0].re,cf0[1].re,cf0[2].re,cf0[3].re);
  // printf("g(x)=%+.5e%+.5e*sqrt(x-a)%+.5e*(x-a)%+.5e*(x-a)**1.5\n",cg0[0].re,cg0[1].re,cg0[2].re,cg0[3].re);
  // printf("%+.5e%+.5e*sqrt(%.5e-x)%+.5e*(%.5e-x)%+.5e*sqrt(%.5e-x)*(%.5e-x)\n",cf[0].im,cf[1].im,METAP_M_MPION_SQUARED,cf[2].im,METAP_M_MPION_SQUARED,cf[3].im,METAP_M_MPION_SQUARED,METAP_M_MPION_SQUARED);
  // printf("%+.5e%+.5e*sqrt(x-%.5e)%+.5e*(x-%.5e)%+.5e*sqrt(x-%.5e)*(x-%.5e)\n",cg[0].im,cg[1].im,METAP_M_MPION_SQUARED,cg[2].im,METAP_M_MPION_SQUARED,cg[3].im,METAP_M_MPION_SQUARED,METAP_M_MPION_SQUARED);
  // printf("a=%.10e\n",METAP_M_MPION_SQUARED);

  si = s[0];
  s_am = s[N[0] + N[1] - 1];
  s_ap = s[N[0] + N[1]];
  sf = s[N[0] + N[1] + N[2] + N[3] - 1];

  s[0] = s0;
  s[N[0] + N[1] - 1] = METAP_M_MPION_SQUARED;
  s[N[0] + N[1]] = METAP_M_MPION_SQUARED;
  s[N[0] + N[1] + N[2] + N[3] - 1] = L2;

  gsl_spline_init(M0[0].re, s, M0_re, N[0] + N[1]);
  gsl_spline_init(M0[0].im, s, M0_im, N[0] + N[1]);

  gsl_spline_init(M0[1].re, s + N[0] + N[1], M0_re + N[0] + N[1], N[2] + N[3]);
  gsl_spline_init(M0[1].im, s + N[0] + N[1], M0_im + N[0] + N[1], N[2] + N[3]);

  gsl_spline_init(M1[0].re, s, M1_re, N[0] + N[1]);
  gsl_spline_init(M1[0].im, s, M1_im, N[0] + N[1]);

  gsl_spline_init(M1[1].re, s + N[0] + N[1], M1_re + N[0] + N[1], N[2] + N[3]);
  gsl_spline_init(M1[1].im, s + N[0] + N[1], M1_im + N[0] + N[1], N[2] + N[3]);

  gsl_spline_init(M2[0].re, s, M2_re, N[0] + N[1]);
  gsl_spline_init(M2[0].im, s, M2_im, N[0] + N[1]);

  gsl_spline_init(M2[1].re, s + N[0] + N[1], M2_re + N[0] + N[1], N[2] + N[3]);
  gsl_spline_init(M2[1].im, s + N[0] + N[1], M2_im + N[0] + N[1], N[2] + N[3]);

  s[0] = si;
  s[N[0] + N[1] - 1] = s_am;
  s[N[0] + N[1]] = s_ap;
  s[N[0] + N[1] + N[2] + N[3] - 1] = sf;

  free(M0_re);
  free(M0_im);

  free(M1_re);
  free(M1_im);

  free(M2_re);
  free(M2_im);
}

void build_inhomogenity_integrand_etap_eta2pi(
    complex_spline *M0, complex_spline *M1, complex *cf0, complex *cg0,
    complex *cf1, complex *cg1, complex *M0_avg, complex *M1_avg_s,
    complex *M1_avg_u, complex *omnes0, complex *omnes1, double *s, double *t,
    double s_th, double s_ps, double a_s, double b_s, double t_th, double t_ps,
    double a_t, double b_t, double L2, int *N_s, int *N_t, int Nsin, int n0,
    int n1) {
  int i, j, n, m, N1, N2;
  double si, sf, s_am, s_ap, k, c, *M0_re, *M0_im, *M1_re, *M1_im, s_f[4],
      s_g[4], abs_omnes;
  complex f0[4], g0[4], f1[4], g1[4];

  M0_re =
      (double *)malloc((N_s[0] + N_s[1] + N_s[2] + N_s[3]) * sizeof(double));
  M0_im =
      (double *)malloc((N_s[0] + N_s[1] + N_s[2] + N_s[3]) * sizeof(double));

  M1_re =
      (double *)malloc((N_t[0] + N_t[1] + N_t[2] + N_t[3]) * sizeof(double));
  M1_im =
      (double *)malloc((N_t[0] + N_t[1] + N_t[2] + N_t[3]) * sizeof(double));

  n = 0;
  for (j = 0; j < 4; j++) {
    if (j != 0) {
      n += N_s[j - 1];
    }

    for (i = 0; i < N_s[j]; i++) {
      k = sqrt(fabs((s[i + n] - s_ps) * (s[i + n] - s_th))) *
          sqrt(fabs(s[i + n] - b_s)) / s[i + n];
      abs_omnes = sqrt(omnes0[i + n].re * omnes0[i + n].re +
                       omnes0[i + n].im * omnes0[i + n].im);
      c = sin(delta0(s[i + n])) / (pow(s[i + n], (double)n0) * k * abs_omnes);

      if (j == 0 || j == 1) {
        M0_re[i + n] = 2. * c * M0_avg[i + n].re;
        M0_im[i + n] = 2. * c * M0_avg[i + n].im;
      }

      else if (j == 2) {
        M0_re[i + n] = 2. * c * M0_avg[i + n].im;
        M0_im[i + n] = -2. * c * M0_avg[i + n].re;
      }

      else if (j == 3) {
        M0_re[i + n] = -2. * c * M0_avg[i + n].re;
        M0_im[i + n] = -2. * c * M0_avg[i + n].im;
      }

      if (isnan(M0_re[i + n]) == 1) {
        printf("FATAL ERROR: M0_hat.re=nan, step=%d,%d\n", i, j);
        exit(1);
      }
      if (isnan(M0_im[i + n]) == 1) {
        printf("FATAL ERROR: M0_hat.im=nan, step=%d,%d\n", i, j);
        exit(1);
      }

      // printf("%.10e %.10e %.10e\n",s[i+n],M0_re[i+n],M0_im[i+n]);
    }
  }

  n = 0;
  for (j = 0; j < 4; j++) {
    if (j == 0) {
      n = 0;
    } else {
      n += N_t[j - 1];
    }

    for (i = 0; i < N_t[j]; i++) {
      k = sqrt(fabs((t[i + n] - t_ps) * (t[i + n] - t_th))) *
          sqrt(fabs(t[i + n] - b_t)) / t[i + n];
      abs_omnes = sqrt(omnes1[i + n].re * omnes1[i + n].re +
                       omnes1[i + n].im * omnes1[i + n].im);
      c = sin(delta_etapi(t[i + n])) /
          (pow(t[i + n], (double)n1) * k * abs_omnes);

      if (j == 0 || j == 1) {
        M1_re[i + n] = c * (-M1_avg_s[i + n].re + M1_avg_u[i + n].re);
        M1_im[i + n] = c * (-M1_avg_s[i + n].im + M1_avg_u[i + n].im);
      }

      else if (j == 2) {
        M1_re[i + n] = c * (-M1_avg_s[i + n].im + M1_avg_u[i + n].im);
        M1_im[i + n] = -c * (-M1_avg_s[i + n].re + M1_avg_u[i + n].re);
      }

      else if (j == 3) {
        M1_re[i + n] = -c * (-M1_avg_s[i + n].re + M1_avg_u[i + n].re);
        M1_im[i + n] = -c * (-M1_avg_s[i + n].im + M1_avg_u[i + n].im);
      }

      if (isnan(M1_re[i + n]) == 1) {
        printf("FATAL ERROR: M1_hat.re=nan, step=%d,%d\n", i, j);
        exit(1);
      }
      if (isnan(M1_im[i + n]) == 1) {
        printf("FATAL ERROR: M1_hat.im=nan, step=%d,%d\n", i, j);
        exit(1);
      }
      //            printf("%.10e %.10e %.10e %.10e %.10e %.10e %.10e
      //            %.10e\n",t[i+n],M1_avg_s[i+n].re,M1_avg_s[i+n].im,M1_avg_u[i+n].re,M1_avg_u[i+n].im,abs_omnes,omnes1[i+n].re,omnes1[i+n].im);

      //            printf("%.10e %.10e %.10e\n",t[i+n],M1_re[i+n],M1_im[i+n]);
    }
  }

  m = 0;
  N1 = N_s[0] + N_s[1] - 1;
  N2 = N_s[0] + N_s[1];
  for (i = 0; i < 4; i++) {
    f0[i].re = M0_re[N1 - m];
    f0[i].im = M0_im[N1 - m];

    g0[i].re = M0_re[N2 + i];
    g0[i].im = M0_im[N2 + i];

    s_f[i] = sqrt(fabs(s[N1 - m] - a_s));
    s_g[i] = sqrt(fabs(s[N2 + i] - a_s));
    m++;
    if (i == 0) {
      N1 -= Nsin;
      N2 += Nsin;
    }
  }

  singularity_sqrt(f0, s_f, cf0);
  singularity_sqrt(g0, s_g, cg0);

  m = 0;
  N1 = N_t[0] + N_t[1] - 1;
  N2 = N_t[0] + N_t[1];
  for (i = 0; i < 4; i++) {
    f1[i].re = M1_re[N1 - m];
    f1[i].im = M1_im[N1 - m];

    g1[i].re = M1_re[N2 + i];
    g1[i].im = M1_im[N2 + i];

    s_f[i] = sqrt(fabs(t[N1 - m] - a_t));
    s_g[i] = sqrt(fabs(t[N2 + i] - a_t));
    m++;
    if (i == 0) {
      N1 -= Nsin;
      N2 += Nsin;
    }
  }

  singularity_sqrt(f1, s_f, cf1);
  singularity_sqrt(g1, s_g, cg1);

  // printf("fre(x)=%+.5e%+.5e*(a-x)%+.5e*(a-x)**1.5%+.5e*(a-x)**2.0%+.5e*(a-x)**2.5\n",cf1[0].re,cf1[1].re,cf1[2].re,cf1[3].re,cf1[4].re);
  // printf("fim(x)=%+.5e%+.5e*(a-x)%+.5e*(a-x)**1.5%+.5e*(a-x)**2.0%+.5e*(a-x)**2.5\n",cf1[0].im,cf1[1].im,cf1[2].im,cf1[3].im,cf1[4].im);
  // printf("gre(x)=%+.5e%+.5e*(x-a)%+.5e*(x-a)**1.5%+.5e*(x-a)**2.0%+.5e*(x-a)**2.5\n",cg1[0].re,cg1[1].re,cg1[2].re,cg1[3].re,cg1[4].re);
  // printf("gim(x)=%+.5e%+.5e*(x-a)%+.5e*(x-a)**1.5%+.5e*(x-a)**2.0%+.5e*(x-a)**2.5\n",cg1[0].im,cg1[1].im,cg1[2].im,cg1[3].im,cg1[4].im);
  //
  //    printf("f(x)=%+.5e%+.5e*sqrt(%.5e-x)%+.5e*(%.5e-x)%+.5e*(%.5e-x)**1.5\n",cf0[0].re,cf0[1].re,a_s,cf0[2].re,a_s,cf0[3].re,a_s);
  //    printf("g(x)=%+.5e%+.5e*sqrt(x-%.5e)%+.5e*(x-%.5e)%+.5e*(x-%.5e)**1.5\n",cg0[0].re,cg0[1].re,a_s,cg0[2].re,a_s,cg0[3].re,a_s);
  //
  //    printf("f(x)=%+.5e%+.5e*sqrt(%.5e-x)%+.5e*(%.5e-x)%+.5e*(%.5e-x)**1.5\n",cf1[0].re,cf1[1].re,a_t,cf1[2].re,a_t,cf1[3].re,a_t);
  //    printf("g(x)=%+.5e%+.5e*sqrt(x-%.5e)%+.5e*(x-%.5e)%+.5e*(x-%.5e)**1.5\n",cg1[0].re,cg1[1].re,a_t,cg1[2].re,a_t,cg1[3].re,a_t);

  // printf("%+.5e%+.5e*sqrt(%.5e-x)%+.5e*(%.5e-x)%+.5e*sqrt(%.5e-x)*(%.5e-x)\n",cf[0].im,cf[1].im,METAP_M_MPION_SQUARED,cf[2].im,METAP_M_MPION_SQUARED,cf[3].im,METAP_M_MPION_SQUARED,METAP_M_MPION_SQUARED);
  // printf("%+.5e%+.5e*sqrt(x-%.5e)%+.5e*(x-%.5e)%+.5e*sqrt(x-%.5e)*(x-%.5e)\n",cg[0].im,cg[1].im,METAP_M_MPION_SQUARED,cg[2].im,METAP_M_MPION_SQUARED,cg[3].im,METAP_M_MPION_SQUARED,METAP_M_MPION_SQUARED);
  // printf("a=%.10e\n",METAP_M_MPION_SQUARED);

  si = s[0];
  s_am = s[N_s[0] + N_s[1] - 1];
  s_ap = s[N_s[0] + N_s[1]];
  sf = s[N_s[0] + N_s[1] + N_s[2] + N_s[3] - 1];

  s[0] = s_th;
  s[N_s[0] + N_s[1] - 1] = a_s;
  s[N_s[0] + N_s[1]] = a_s;
  s[N_s[0] + N_s[1] + N_s[2] + N_s[3] - 1] = L2;

  gsl_spline_init(M0[0].re, s, M0_re, N_s[0] + N_s[1]);
  gsl_spline_init(M0[0].im, s, M0_im, N_s[0] + N_s[1]);

  gsl_spline_init(M0[1].re, s + N_s[0] + N_s[1], M0_re + N_s[0] + N_s[1],
                  N_s[2] + N_s[3]);
  gsl_spline_init(M0[1].im, s + N_s[0] + N_s[1], M0_im + N_s[0] + N_s[1],
                  N_s[2] + N_s[3]);

  s[0] = si;
  s[N_s[0] + N_s[1] - 1] = s_am;
  s[N_s[0] + N_s[1]] = s_ap;
  s[N_s[0] + N_s[1] + N_s[2] + N_s[3] - 1] = sf;

  si = t[0];
  s_am = t[N_t[0] + N_t[1] - 1];
  s_ap = t[N_t[0] + N_t[1]];
  sf = t[N_t[0] + N_t[1] + N_t[2] + N_t[3] - 1];

  t[0] = t_th;
  t[N_t[0] + N_t[1] - 1] = a_t;
  t[N_t[0] + N_t[1]] = a_t;
  t[N_t[0] + N_t[1] + N_t[2] + N_t[3] - 1] = L2;

  gsl_spline_init(M1[0].re, t, M1_re, N_t[0] + N_t[1]);
  gsl_spline_init(M1[0].im, t, M1_im, N_t[0] + N_t[1]);

  gsl_spline_init(M1[1].re, t + N_t[0] + N_t[1], M1_re + N_t[0] + N_t[1],
                  N_t[2] + N_t[3]);
  gsl_spline_init(M1[1].im, t + N_t[0] + N_t[1], M1_im + N_t[0] + N_t[1],
                  N_t[2] + N_t[3]);

  t[0] = si;
  t[N_t[0] + N_t[1] - 1] = s_am;
  t[N_t[0] + N_t[1]] = s_ap;
  t[N_t[0] + N_t[1] + N_t[2] + N_t[3] - 1] = sf;

  free(M0_re);
  free(M0_im);

  free(M1_re);
  free(M1_im);
}

void build_inhomogenity_integrand_Matrix(
    complex_spline *M0, complex_spline *M1, complex_spline *M2, complex *cf0,
    complex *cg0, complex *cf1, complex *cg1, complex *cf2, complex *cg2,
    complex *M0_tilde, complex *M1_tilde, complex *M2_tilde, double *sin_delta0,
    double *sin_delta1, double *sin_delta2, double *abs_omnes0,
    double *abs_omnes1, double *abs_omnes2, double *s, double s0, double L2,
    int *N, int n0, int n1, int n2) {
  int i, n, m, N1, N2, low[5], high[5];
  double si, sf, s_am, a, b, s_ap, ds, k, c0, c1, c2, *M0_re, *M0_im, *M1_re,
      *M1_im, *M2_re, *M2_im, s_f[5], s_g[5];
  complex f0[4], g0[4], f1[5], g1[5], f2[4], g2[4];

  a = METAP_M_MPION_SQUARED;
  b = METAP_P_MPION_SQUARED;

  n = N[0] + N[1];

  M0_re = (double *)malloc(n * sizeof(double));
  M0_im = (double *)malloc(n * sizeof(double));

  M1_re = (double *)malloc(n * sizeof(double));
  M1_im = (double *)malloc(n * sizeof(double));

  M2_re = (double *)malloc(n * sizeof(double));
  M2_im = (double *)malloc(n * sizeof(double));

  ds = .25;
  if (ds < 0.5 * (s[N[0] - 1] - s[N[0] - 2] + s[N[0] + 1] - s[N[0]])) {
    printf("FATAL ERROR: In function 'build_inhomogenity_integrand_Matrix', s "
           "stepsize is to large for calculating fit-functions at "
           "pseudo-threshold! ds=%.3e>%.3e\n\n",
           0.5 * (s[N[0] - 1] - s[N[0] - 2] + s[N[0] + 1] - s[N[0]]), ds);
    exit(1);
  }

  for (i = 0; i < n; i++) {

    k = sqrt(1. - 4. / s[i]) * sqrt(fabs(s[i] - b));

    c0 = sin_delta0[i] / (pow(s[i], (double)n0) * k * abs_omnes0[i]);
    c1 = sin_delta1[i] / (pow(s[i], (double)n1) * k * k * k * abs_omnes1[i]);
    c2 = sin_delta2[i] / (pow(s[i], (double)n2) * k * abs_omnes2[i]);

    if (s[i] < a) {
      M0_re[i] = c0 * M0_tilde[i].re;
      M0_im[i] = c0 * M0_tilde[i].im;
      M1_re[i] = c1 * M1_tilde[i].re;
      M1_im[i] = c1 * M1_tilde[i].im;
      M2_re[i] = c2 * M2_tilde[i].re;
      M2_im[i] = c2 * M2_tilde[i].im;
    } else if (s[i] > a && s[i] < b) {
      M0_re[i] = c0 * M0_tilde[i].im;
      M0_im[i] = -c0 * M0_tilde[i].re;
      M1_re[i] = -c1 * M1_tilde[i].im;
      M1_im[i] = c1 * M1_tilde[i].re;
      M2_re[i] = c2 * M2_tilde[i].im;
      M2_im[i] = -c2 * M2_tilde[i].re;
    } else {
      M0_re[i] = -c0 * M0_tilde[i].re;
      M0_im[i] = -c0 * M0_tilde[i].im;
      M1_re[i] = -c1 * M1_tilde[i].re;
      M1_im[i] = -c1 * M1_tilde[i].im;
      M2_re[i] = -c2 * M2_tilde[i].re;
      M2_im[i] = -c2 * M2_tilde[i].im;
    }

    if (s[i] < a - 8. * ds) {
      low[4] = i;
    }
    if (s[i] < a + 8. * ds) {
      high[4] = i;
    }
    if (s[i] < a - 4. * ds) {
      low[3] = i;
    }
    if (s[i] < a + 4. * ds) {
      high[3] = i;
    }
    if (s[i] < a - 2. * ds) {
      low[2] = i;
    }
    if (s[i] < a + 2. * ds) {
      high[2] = i;
    }
    if (s[i] < a - ds) {
      low[1] = i;
    }
    if (s[i] < a + ds) {
      high[1] = i;
    }

    // printf("%.3e %.3e %.3e %.3e %.3e %.3e
    // %.3e\n",s[i],M0_re[i],M0_im[i],M1_re[i],M1_im[i],M2_re[i],M2_im[i]);

    if (isnan(M0_re[i]) == 1) {
      printf("FATAL ERROR: M0_hat.re=nan, step=%d\n", i);
      printf("%.3e %.3e %.3e %.3e %.3e\n", s[i], M0_tilde[i].re, sin_delta0[i],
             abs_omnes0[i], c0);
      exit(1);
    }
    if (isnan(M0_im[i]) == 1) {
      printf("FATAL ERROR: M0_hat.im=nan, step=%d\n", i);
      exit(1);
    }
    if (isnan(M1_re[i]) == 1) {
      printf("FATAL ERROR: M1_hat.re=nan, step=%d\n", i);
      exit(1);
    }
    if (isnan(M1_im[i]) == 1) {
      printf("FATAL ERROR: M1_hat.im=nan, step=%d\n", i);
      exit(1);
    }
    if (isnan(M2_re[i]) == 1) {
      printf("FATAL ERROR: M2_hat.re=nan, step=%d\n", i);
      exit(1);
    }
    if (isnan(M2_im[i]) == 1) {
      printf("FATAL ERROR: M2_hat.im=nan, step=%d\n", i);
      exit(1);
    }
  }

  low[0] = N[0] - 1;
  high[0] = N[0];

  for (i = 0; i < 4; i++) {
    f0[i].re = M0_re[low[i]];
    f0[i].im = M0_im[low[i]];

    f1[i].re = M1_re[low[i]];
    f1[i].im = M1_im[low[i]];

    f2[i].re = M2_re[low[i]];
    f2[i].im = M2_im[low[i]];

    g0[i].re = M0_re[high[i]];
    g0[i].im = M0_im[high[i]];

    g1[i].re = M1_re[high[i]];
    g1[i].im = M1_im[high[i]];

    g2[i].re = M2_re[high[i]];
    g2[i].im = M2_im[high[i]];

    s_f[i] = sqrt(fabs(s[low[i]] - a));
    s_g[i] = sqrt(fabs(s[high[i]] - a));
  }
  f1[i].re = M1_re[low[i]];
  f1[i].im = M1_im[low[i]];

  g1[i].re = M1_re[high[i]];
  g1[i].im = M1_im[high[i]];

  s_f[i] = sqrt(fabs(s[low[i]] - a));
  s_g[i] = sqrt(fabs(s[high[i]] - a));

  singularity_sqrt(f0, s_f, cf0);
  singularity_sqrt(g0, s_g, cg0);

  singularity_sqrt_cubed(f1, s_f, cf1);
  singularity_sqrt_cubed(g1, s_g, cg1);

  singularity_sqrt(f2, s_f, cf2);
  singularity_sqrt(g2, s_g, cg2);

  si = s[0];
  s_am = s[N[0] - 1];
  s_ap = s[N[0]];
  sf = s[N[0] + N[1] - 1];

  s[0] = s0;
  s[N[0] - 1] = a;
  s[N[0]] = a;
  s[N[0] + N[1] - 1] = L2;

  gsl_spline_init(M0[0].re, s, M0_re, N[0]);
  gsl_spline_init(M0[0].im, s, M0_im, N[0]);

  gsl_spline_init(M0[1].re, s + N[0], M0_re + N[0], N[1]);
  gsl_spline_init(M0[1].im, s + N[0], M0_im + N[0], N[1]);

  gsl_spline_init(M1[0].re, s, M1_re, N[0]);
  gsl_spline_init(M1[0].im, s, M1_im, N[0]);

  gsl_spline_init(M1[1].re, s + N[0], M1_re + N[0], N[1]);
  gsl_spline_init(M1[1].im, s + N[0], M1_im + N[0], N[1]);

  gsl_spline_init(M2[0].re, s, M2_re, N[0]);
  gsl_spline_init(M2[0].im, s, M2_im, N[0]);

  gsl_spline_init(M2[1].re, s + N[0], M2_re + N[0], N[1]);
  gsl_spline_init(M2[1].im, s + N[0], M2_im + N[0], N[1]);

  s[0] = si;
  s[N[0] - 1] = s_am;
  s[N[0]] = s_ap;
  s[N[0] + N[1] - 1] = sf;

  free(M0_re);
  free(M0_im);

  free(M1_re);
  free(M1_im);

  free(M2_re);
  free(M2_im);
}

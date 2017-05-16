#include "Amplitude.h"
#include "Basic.h"

void subtraction_constant(char *sub_const, int *n0, int *n1, int *n2) {
  if (str_eq_str(sub_const, "a0") == 1) {
    *n0 = 0;
    *n1 = -1;
    *n2 = -1;
  } else if (str_eq_str(sub_const, "b0") == 1) {
    *n0 = 1;
    *n1 = -1;
    *n2 = -1;
  } else if (str_eq_str(sub_const, "c0") == 1) {
    *n0 = 2;
    *n1 = -1;
    *n2 = -1;
  } else if (str_eq_str(sub_const, "d0") == 1) {
    *n0 = 3;
    *n1 = -1;
    *n2 = -1;
  } else if (str_eq_str(sub_const, "a1") == 1) {
    *n0 = -1;
    *n1 = 0;
    *n2 = -1;
  } else if (str_eq_str(sub_const, "b1") == 1) {
    *n0 = -1;
    *n1 = 1;
    *n2 = -1;
  } else if (str_eq_str(sub_const, "c1") == 1) {
    *n0 = -1;
    *n1 = 2;
    *n2 = -1;
  } else if (str_eq_str(sub_const, "a2") == 1) {
    *n0 = -1;
    *n1 = -1;
    *n2 = 0;
  } else if (str_eq_str(sub_const, "b2") == 1) {
    *n0 = -1;
    *n1 = -1;
    *n2 = 1;
  } else if (str_eq_str(sub_const, "etapi") == 1) {
    *n0 = -1;
    *n1 = -1;
    *n2 = -1;
  } else if (str_eq_str(sub_const, "test") == 1) {
    *n0 = 0;
    *n1 = 1;
    *n2 = -1;
  } else {
    printf("FATAL ERROR: Bad input, unknown subtraction constant: %s\n",
           sub_const);
    exit(1);
  }
}

void build_inhomogenity_init(complex **M_inhom, int *N) {
  int i, j;

  for (i = 0; i < 7; i++) {
    for (j = 0; j < N[i]; j++) {
      M_inhom[i][j].re = 0.0;
      M_inhom[i][j].im = 0.0;
    }
  }
}

void amplitude_cv_plus(complex_spline M, complex *omnes, complex *M_inhom,
                       double *s, double s_th, double L2, int N, int n) {
  int i;
  double si, sf;
  complex m, O;
  double *M_re = (double *)malloc(N * sizeof(double));
  double *M_im = (double *)malloc(N * sizeof(double));

  for (i = 0; i < N; i++) {
    O.re = omnes[i].re;
    O.im = omnes[i].im;
    m.re = M_inhom[i].re;
    m.im = M_inhom[i].im;

    if (isnan(m.re) == 1) {
      printf("FATAL ERROR: In function amplitude_cv_plus: m.re=nan, step=%d\n",
             i);
      exit(1);
    } else if (isnan(m.im) == 1) {
      printf("FATAL ERROR: In function amplitude_cv_plus: m.im=nan, step=%d\n",
             i);
      exit(1);
    }

    if (isnan(O.re) == 1) {
      printf("FATAL ERROR: In function amplitude_cv_plus: O.re=nan, step=%d\n",
             i);
      exit(1);
    } else if (isnan(O.im) == 1) {
      printf("FATAL ERROR: In function amplitude_cv_plus: O.im=nan, step=%d\n",
             i);
      exit(1);
    }

    if (n == -1) {
      M_re[i] = O.re * m.re - O.im * m.im;
      M_im[i] = O.re * m.im + O.im * m.re;
    }

    else if (n >= 0) {
      M_re[i] = O.re * (pow(s[i], n) + m.re) - O.im * m.im;
      M_im[i] = O.re * m.im + O.im * (pow(s[i], n) + m.re);
    }
  }

  si = s[0];
  sf = s[N - 1];
  s[0] = s_th;
  s[N - 1] = L2;
  gsl_spline_init(M.re, s, M_re, N);
  gsl_spline_init(M.im, s, M_im, N);
  s[0] = si;
  s[N - 1] = sf;

  free(M_re);
  free(M_im);
}

void amplitude_cv_minus(complex_spline M, complex *omnes, complex *M_inhom,
                        double *s, double s_th, int N, int n) {
  int i;
  double si;
  complex m, O;
  double *M_re = (double *)malloc(N * sizeof(double));
  double *M_im = (double *)malloc(N * sizeof(double));

  for (i = 0; i < N; i++) {
    O.re = omnes[i].re;
    O.im = omnes[i].im;
    m.re = M_inhom[i].re;
    m.im = M_inhom[i].im;

    if (isnan(m.re) == 1) {
      printf("FATAL ERROR: In function amplitude_cv_minus: m.re=nan, step=%d\n",
             i);
      exit(1);
    } else if (isnan(m.im) == 1) {
      printf("FATAL ERROR: In function amplitude_cv_minus: m.im=nan, step=%d\n",
             i);
      exit(1);
    }

    if (isnan(O.re) == 1) {
      printf("FATAL ERROR: In function amplitude_cv_minus: O.re=nan, step=%d\n",
             i);
      exit(1);
    } else if (isnan(O.im) == 1) {
      printf("FATAL ERROR: In function amplitude_cv_minus: O.im=nan, step=%d\n",
             i);
      exit(1);
    }

    if (n == -1) {
      M_re[i] = O.re * m.re - O.im * m.im;
      M_im[i] = O.re * m.im + O.im * m.re;
    }

    else if (n >= 0) {
      M_re[i] = O.re * (pow(s[i], n) + m.re) - O.im * m.im;
      M_im[i] = O.re * m.im + O.im * (pow(s[i], n) + m.re);
    }
  }

  si = s[0];
  s[0] = s_th;
  gsl_spline_init(M.re, s, M_re, N);
  gsl_spline_init(M.im, s, M_im, N);
  s[0] = si;

  free(M_re);
  free(M_im);
}

void amplitude_below_cut(complex_spline M, complex *omnes, complex *M_inhom,
                         double *s, double s_th, int N, int n) {
  int i;
  double sf;
  complex m, O;
  double *M_re = (double *)malloc(N * sizeof(double));
  double *M_im = (double *)malloc(N * sizeof(double));

  for (i = 0; i < N; i++) {
    O.re = omnes[i].re;
    O.im = omnes[i].im;
    m.re = M_inhom[i].re;
    m.im = M_inhom[i].im;

    if (isnan(m.re) == 1) {
      printf(
          "FATAL ERROR: In function amplitude_below_cut: m.re=nan, step=%d\n",
          i);
      exit(1);
    } else if (isnan(m.im) == 1) {
      printf(
          "FATAL ERROR: In function amplitude_below_cut: m.im=nan, step=%d\n",
          i);
      exit(1);
    }

    if (isnan(O.re) == 1) {
      printf(
          "FATAL ERROR: In function amplitude_below_cut: O.re=nan, step=%d\n",
          i);
      exit(1);
    } else if (isnan(O.im) == 1) {
      printf(
          "FATAL ERROR: In function amplitude_below_cut: O.im=nan, step=%d\n",
          i);
      exit(1);
    }

    if (n == -1) {
      M_re[i] = O.re * m.re - O.im * m.im;
      M_im[i] = O.re * m.im + O.im * m.re;
    }

    else if (n >= 0) {
      M_re[i] = O.re * (pow(s[i], n) + m.re) - O.im * m.im;
      M_im[i] = O.re * m.im + O.im * (pow(s[i], n) + m.re);
    }
  }

  sf = s[N - 1];
  s[N - 1] = s_th;
  gsl_spline_init(M.re, s, M_re, N);
  gsl_spline_init(M.im, s, M_im, N);
  s[N - 1] = sf;

  free(M_re);
  free(M_im);
}

void amplitude_complex(complex_spline M, complex *omnes, complex *M_inhom,
                       complex *s, double *y, char plusminus, int N, int n) {
  int i;
  double r, phi;
  complex m, O;
  double *M_re = (double *)malloc(N * sizeof(double));
  double *M_im = (double *)malloc(N * sizeof(double));

  for (i = 0; i < N; i++) {
    O.re = omnes[i].re;
    O.im = omnes[i].im;
    m.re = M_inhom[i].re;
    m.im = M_inhom[i].im;

    switch (plusminus) {
    case '+':
      break;

    case '-':
      O.im *= -1.;
      s[i].im *= -1.;
      break;

    default:
      printf("FATAL ERROR: In function inhomogenity_sqrt_complex, plusminus "
             "has to be '+' or '-', '%c' is unknown!",
             plusminus);
      exit(1);
      break;
    }

    if (isnan(m.re) == 1) {
      printf("FATAL ERROR: In function amplitude_complex: m.re=nan, step=%d\n",
             i);
      exit(1);
    } else if (isnan(m.im) == 1) {
      printf("FATAL ERROR: In function amplitude_complex: m.im=nan, step=%d\n",
             i);
      exit(1);
    }

    if (isnan(O.re) == 1) {
      printf("FATAL ERROR: In function amplitude_complex: O.re=nan, step=%d\n",
             i);
      exit(1);
    } else if (isnan(O.im) == 1) {
      printf("FATAL ERROR: In function amplitude_complex: O.im=nan, step=%d\n",
             i);
      exit(1);
    }

    if (n == -1) {
      M_re[i] = O.re * m.re - O.im * m.im;
      M_im[i] = O.re * m.im + O.im * m.re;
    }

    else if (n >= 0) {
      r = pow(sqrt(s[i].re * s[i].re + s[i].im * s[i].im), n);
      phi = n * atan2(s[i].im, s[i].re);

      M_re[i] = O.re * (r * cos(phi) + m.re) - O.im * (r * sin(phi) + m.im);
      M_im[i] = O.re * (r * sin(phi) + m.im) + O.im * (r * cos(phi) + m.re);
    }

    switch (plusminus) {
    case '+':
      break;

    case '-':
      s[i].im *= -1.;
      break;

    default:
      printf("FATAL ERROR: In function inhomogenity_sqrt_complex, plusminus "
             "has to be '+' or '-', '%c' is unknown!",
             plusminus);
      exit(1);
      break;
    }
  }

  gsl_spline_init(M.re, y, M_re, N);
  gsl_spline_init(M.im, y, M_im, N);

  free(M_re);
  free(M_im);
}

void build_amplitude(complex_spline *M, complex **omnes, complex **M_inhom,
                     double *s_cv_plus, double *s_cv_minus, double *s_below_cut,
                     complex **s_cmp, double **y, double s0, double L2, int *N,
                     int n) {

  amplitude_cv_plus(M[0], omnes[0], M_inhom[0], s_cv_plus, s0, L2, N[0], n);
  amplitude_cv_minus(M[1], omnes[1], M_inhom[1], s_cv_minus, s0, N[1], n);
  amplitude_below_cut(M[2], omnes[2], M_inhom[2], s_below_cut, s0, N[2], n);
  amplitude_complex(M[3], omnes[3], M_inhom[3], s_cmp[0], y[0], '-', N[3], n);
  amplitude_complex(M[4], omnes[4], M_inhom[4], s_cmp[1], y[1], '-', N[3], n);
  amplitude_complex(M[5], omnes[3], M_inhom[5], s_cmp[0], y[0], '+', N[3], n);
  amplitude_complex(M[6], omnes[4], M_inhom[6], s_cmp[1], y[1], '+', N[3], n);
}

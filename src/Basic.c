#include "Basic.h"

// MKAON^2 in MPION^2
const double MKAON_SQUARED = MKAON * MKAON / (MPION * MPION);

// META^2 in MPION^2
const double META_SQUARED = META * META / (MPION * MPION);

// MRHO^2 in MPION^2
const double MRHO_SQUARED = MRHO * MRHO / (MPION * MPION);

// METAP^2 in MPION^2
const double METAP_SQUARED = METAP * METAP / (MPION * MPION);

//(METAP-MPION)^2 in MPION^2
const double METAP_M_MPION_SQUARED =
    (METAP / MPION - 1.0) * (METAP / MPION - 1.0);

//(METAP+MPION)^2 in MPION^2
const double METAP_P_MPION_SQUARED =
    (METAP / MPION + 1.0) * (METAP / MPION + 1.0);

//(META+MPION)^2 in MPION^2
const double s_etapi = (META / MPION + 1.) * (META / MPION + 1.);

//(META-MPION)^2 in MPION^2
const double s_etapi_pseudo = (META / MPION - 1.) * (META / MPION - 1.);

// real part of complex variable
double cmp_re(complex cmp) { return cmp.re; }

// imaginary part of complex variable
double cmp_im(complex cmp) { return cmp.im; }

// Phase space for etapi-disc
double lambda_etapi(double s) {
  if (s <= s_etapi) {
    return 0.;
  } else {
    return sqrt((s_etapi_pseudo - s) * (s_etapi - s)) / s;
  }
}

// string comparison
int str_eq_str(char *str1, char *str2) {
  int i, n;

  n = strlen(str2);
  for (i = 0; str1[i] == *(str2 + i) && i < n; i++) {
  }

  if (i == n) {
    return 1;
  } else {
    return 0;
  }
}

// 1/x for complex values of x
complex c_inverse(complex x) {
  double temp;
  complex result;

  temp = 1.0 / (x.re * x.re + x.im * x.im);
  result.re = x.re * temp;
  result.im = -x.im * temp;

  return result;
}

// returns radius and phase of complex variable x
void polar_r_phi(complex x, double *r, double *phi) {
  *r = sqrt(x.re * x.re + x.im * x.im);
  *phi = atan2(x.im, x.re);
}

// sqrt(x) for complex values of x
complex c_sqrt(complex x) {
  double r, phi, cphi, sphi;
  complex result;

  polar_r_phi(x, &r, &phi);

  r = sqrt(r);
  phi = 0.5 * phi;

  cphi = cos(phi);
  sphi = sin(phi);

  result.re = r * cphi;
  result.im = r * sphi;

  return result;
}

// log(x) for complex values of x
complex c_log(complex x) {
  double r, phi;
  complex result;

  polar_r_phi(x, &r, &phi);

  result.re = log(r);
  result.im = phi;

  return result;
}

// 2*atanh(x) for complex values of x
complex c_two_atanh(complex x) {
  double temp;
  complex result;

  temp = 1.0 - x.re;
  temp = 1.0 / (temp * temp + x.im * x.im);

  result.re = (1.0 - x.re * x.re - x.im * x.im) * temp;
  result.im = 2.0 * x.im * temp;

  result = c_log(result);

  return result;
}

// 2*i*atan(x) for complex values of x
complex c_two_i_atan(complex x) {
  double temp;
  complex result;

  temp = 1.0 + x.im;
  temp = 1.0 / (temp * temp + x.re * x.re);

  result.re = (1.0 - x.re * x.re - x.im * x.im) * temp;
  result.im = 2.0 * x.re * temp;

  result = c_log(result);

  return result;
}

complex kappa(double s) {
  complex result;
  double temp;

  temp = (1.0 - 4.0 / s) * (METAP_M_MPION_SQUARED - s) *
         (METAP_P_MPION_SQUARED - s);

  if (temp > 0.0) {
    result.re = sqrt(temp);
    result.im = 0.0;
    if (s > METAP_P_MPION_SQUARED) {
      result.re *= -1.0;
    }
  }

  else {
    result.re = 0.0;
    result.im = sqrt(-temp);
  }

  return result;
}

void phase_shifts(gsl_spline *delta0, gsl_spline *delta1, gsl_spline *delta2,
                  double *sin_delta0, double *sin_delta1, double *sin_delta2,
                  double *s, int N) {
  int i;

  gsl_interp_accel *acc = gsl_interp_accel_alloc();
  for (i = 0; i < N; i++) {
    sin_delta0[i] = sin(gsl_spline_eval(delta0, s[i], acc));
    sin_delta1[i] = sin(gsl_spline_eval(delta1, s[i], acc));
    sin_delta2[i] = sin(gsl_spline_eval(delta2, s[i], acc));
  }
  gsl_interp_accel_free(acc);
}

// angle phi of the complex path for given s
double integration_III_phi_s(double s) {
  double result;

  if (s == METAP_M_MPION_SQUARED) {
    result = 0.0;
  } else if (s == METAP_P_MPION_SQUARED) {
    result = M_PI;
  } else {
    result =
        acos(0.5 * sqrt(s) * (METAP_SQUARED + 3.0 - s) / (METAP_SQUARED - 1.0));
  }

  return result;
}

// radius of the complex path for given s
double integration_III_radius_s(double s) {
  return sqrt(METAP_M_MPION_SQUARED * METAP_P_MPION_SQUARED / s);
}

// radius of the complex path for given phi
double integration_III_radius_phi(double phi) {
  double result, z, eta, eps;

  eps = 1.0e-05;

  z = cos(phi);

  if (z == -1.0) {
    result = sqrt(METAP_M_MPION_SQUARED);
  } else if (z > -1.0 && z < -eps) {
    eta = 1.0 / 3.0 * (acos(54.0 * pow(METAP_SQUARED - 1.0, 2.0) /
                                pow(METAP_SQUARED + 3.0, 3.0) * pow(z, 2.0) -
                            1.0));
    result = (METAP_SQUARED + 3.0) / (6.0 * z) * (1.0 - 2.0 * cos(eta));
  } else if (z >= -eps && z <= eps) {
    result = (METAP_SQUARED - 1.0) / sqrt(METAP_SQUARED + 3.0);
  } else if (z > eps && z < 1.0) {
    eta = 1.0 / 3.0 * (acos(54.0 * pow(METAP_SQUARED - 1.0, 2.0) /
                                pow(METAP_SQUARED + 3.0, 3.0) * pow(z, 2.0) -
                            1.0) +
                       M_PI);
    result = (METAP_SQUARED + 3.0) / (6.0 * z) * (1.0 + 2.0 * cos(eta));
  } else if (z == 1.0) {
    result = sqrt(METAP_P_MPION_SQUARED);
  }

  return result;
}

// derivative of the radius with respect to phi for given phi
double integration_III_dradius_phi(double phi) {
  double result, z, eta, gamma, eps;

  eps = 1.0e-05;

  z = cos(phi);

  if (z >= -1.0 && z < -eps) {
    gamma = 54.0 * pow(METAP_SQUARED - 1.0, 2.0) /
                pow(METAP_SQUARED + 3.0, 3.0) * pow(z, 2.0) -
            1.0;
    eta = 1.0 / 3.0 * acos(gamma);
    result =
        -(METAP_SQUARED + 3.0) / (6.0 * pow(z, 2.0)) * (1.0 - 2.0 * cos(eta)) -
        12.0 * pow(METAP_SQUARED - 1.0, 2.0) * sin(eta) /
            (pow(METAP_SQUARED + 3.0, 2.0) * sqrt(1.0 - pow(gamma, 2.0)));
  } else if (z >= -eps && z <= eps) {
    result = pow(METAP_SQUARED - 1.0, 2.0) / pow(METAP_SQUARED + 3.0, 2.0);
  } else if (z > eps && z <= 1.0) {
    gamma = 54.0 * pow(METAP_SQUARED - 1.0, 2.0) /
                pow(METAP_SQUARED + 3.0, 3.0) * pow(z, 2.0) -
            1.0;
    eta = 1.0 / 3.0 * (acos(gamma) + M_PI);
    result =
        -(METAP_SQUARED + 3.0) / (6.0 * pow(z, 2.0)) * (1.0 + 2.0 * cos(eta)) +
        12.0 * pow(METAP_SQUARED - 1.0, 2.0) * sin(eta) /
            (pow(METAP_SQUARED + 3.0, 2.0) * sqrt(1.0 - pow(gamma, 2.0)));
  }

  return -result * sin(phi);
}

// path gamma for given phi
complex integration_III_gamma_phi(double phi) {
  double r;
  complex result;

  r = integration_III_radius_phi(phi);

  result.re = r * cos(phi);
  result.im = r * sin(phi);

  return result;
}

// derivative of the path gamma with respect to phi for given phi
complex integration_III_dgamma_phi(double phi) {
  double x, y, r, dr;
  complex result;

  x = sin(phi);
  y = cos(phi);
  r = integration_III_radius_phi(phi);
  dr = integration_III_dradius_phi(phi);

  result.re = -r * x + dr * y;
  result.im = r * y + dr * x;

  return result;
}

double integration_III_circumference() {
  int M = 1500;
  double result, error, eps1, phi;
  eps1 = 1.0e-07;
  double f(double phi, void *params) {
    double r = integration_III_radius_phi(phi);
    double dr = integration_III_dradius_phi(phi);
    double f = sqrt(pow(r, 2.0) + pow(dr, 2.0));
    return f;
  }
  gsl_integration_workspace *w = gsl_integration_workspace_alloc(M);

  gsl_function F;
  F.params = NULL;

  F.function = &f;
  gsl_integration_qags(&F, 0.0 + eps1, 2.0 * M_PI - eps1, eps1, eps1, M, w,
                       &result, &error);
  gsl_integration_workspace_free(w);

  return result;
}

// New integration path III
// calc new integration variable ya or yb from s with given a=(M-m)**2 and
// b=(M+m)**2 in range s in [a,b]
double s_to_y(double s, double a, double b) {
  double a_b_avg = 0.5 * (a + b);

  if (s >= a && s <= a_b_avg) {
    return sqrt(0.5 * (sqrt(s) + sqrt(a))) - sqrt(0.5 * (sqrt(s) - sqrt(a)));
  } else if (s > a_b_avg && s <= b) {
    return sqrt(0.5 * (sqrt(b) + sqrt(s))) - sqrt(0.5 * (sqrt(b) - sqrt(s)));
  } else {
    printf("FATAL ERROR: In function s_to_y, s=%.3e not in range [%.3e,%.3e]!",
           s, a, b);
    exit(1);
  }
}

// calc s from ya or yb with given a=(M-m)**2 and b=(M+m)**2 in range s in [a,b]
double y_to_s(double y, double a, double b, char region) {
  double s, y2, y4;

  switch (region) {
  case 'a':
    y4 = y * y * y * y;
    s = pow(y4 + a, 2.) / (4. * y4);
    break;

  case 'b':
    y2 = y * y;
    s = y2 * (2. * sqrt(b) - y2);
    break;

  default:
    printf("FATAL ERROR: In function y_to_s, unknown region '%c'!", region);
    exit(1);
    break;
  }

  return s;
}

// integration path gamma for s>=d==s_th
complex gamma_path(double s, double a, double b, double c, double d,
                   char plusminus, char kappa_pm) {
  complex gamma, kappa;

  switch (plusminus) {
  case '+':
    gamma.re = 0.5 * (0.5 * (a + b + c + d) - s + sqrt(a * b * c * d) / s);
    break;

  case '0':
    gamma.re = 0.5 * (0.5 * (a + b + c + d) - s);
    break;

  case '-':
    gamma.re = 0.5 * (0.5 * (a + b + c + d) - s - sqrt(a * b * c * d) / s);
    break;

  default:
    printf("FATAL ERROR: In function gamma_path, plusminus has to be '+' or "
           "'0' or '-', '%c' is unknown!",
           plusminus);
    exit(1);
    break;
  }

  if (s == d || s == a || s == b) {
    kappa.re = 0.;
    kappa.im = 0.;
  }

  else if (s > d && s < a) {
    kappa.re = 0.5 * (sqrt((s - c) * (s - d)) * sqrt((s - a) * (s - b))) / s;
    kappa.im = 0.;
  }

  else if (s > a && s < b) {
    kappa.re = 0.;
    kappa.im = 0.5 * (sqrt((s - c) * (s - d)) * sqrt((s - a) * (b - s))) / s;
  }

  else if (s > b) {
    kappa.re = -0.5 * (sqrt((s - c) * (s - d)) * sqrt((s - a) * (s - b))) / s;
    kappa.im = 0.;
  }

  else {
    printf("FATAL ERROR: In function gamma_path, value of s is not covered by "
           "definition range!\n");
    exit(1);
  }

  gamma.im = 0.;

  switch (kappa_pm) {
  case '+':
    gamma.re += kappa.re;
    gamma.im += kappa.im;
    break;

  case '-':
    gamma.re -= kappa.re;
    gamma.im -= kappa.im;
    break;

  default:
    printf("FATAL ERROR: In function gamma_path, kappa_pm has to be '+' or "
           "'-', '%c' is unknown!",
           kappa_pm);
    exit(1);
    break;
  }

  return gamma;
}

// integration path gamma in terms of intrgarion variable y, a=(m1-m2)**2,
// b=(m1+m2)**2, c=(m3-m4)**2, d=(m3+m4)**2
// m1,...,m4 particle masses
complex gamma_III(double y, double a, double b, double c, double d, char region,
                  char plusminus, char kappa_pm) {
  complex gamma;
  double s;

  s = y_to_s(y, a, b, region);

  //    switch (plusminus) {
  //        case '+':
  //            gamma.re = 0.5*(0.5*(a+b+c+d)-s+sqrt(a*b*c*d)/s);
  //            break;
  //
  //        case '0':
  //            gamma.re = 0.5*(0.5*(a+b+c+d)-s);
  //            break;
  //
  //        case '-':
  //            gamma.re = 0.5*(0.5*(a+b+c+d)-s-sqrt(a*b*c*d)/s);
  //            break;
  //
  //        default:
  //            printf("FATAL ERROR: In function gamma_III, plusminus has to be
  //            '+' or '0' or '-', '%c' is unknown!",plusminus);
  //            exit(1);
  //            break;
  //    }
  //
  //    gamma.im = 0.5*(sqrt((s-c)*(s-d))*sqrt((a-s)*(s-b)))/s;

  return gamma_path(s, a, b, c, d, plusminus, kappa_pm);
}

// derivative of the integration path gamma in terms of intrgarion variable y,
// a=(m1-m2)**2, b=(m1+m2)**2, c=(m3-m4)**2, d=(m3+m4)**2
// m1,...,m4 particle masses
complex dgamma_III(double y, double a, double b, double c, double d,
                   char region, char plusminus, char kappa_pm) {
  complex dgamma;
  double s, dsdy;

  s = y_to_s(y, a, b, region);

  switch (plusminus) {
  case '+':
    dgamma.re = -0.5 * (1. + sqrt(a * b * c * d) / (s * s));
    break;

  case '0':
    dgamma.re = -0.5;
    break;

  case '-':
    dgamma.re = -0.5 * (1. - sqrt(a * b * c * d) / (s * s));
    break;

  default:
    printf("FATAL ERROR: In function dgamma_III, plusminus has to be '+' or "
           "'0' or '-', '%c' is unknown!",
           plusminus);
    exit(1);
    break;
  }

  dgamma.im = (2. * a * b * c * d - (b * c * (a + d) + a * d * (b + c)) * s +
               (a + b + c + d) * s * s * s - 2. * s * s * s * s) /
              (4. * s * s * sqrt((s - a) * (b - s)) * sqrt((s - c) * (s - d)));

  switch (region) {
  case 'a':
    dsdy = y * y * y - a * a / pow(y, 5.);

    if (y >= pow(a, 1. / 4.) - 1.e-3) {
      dgamma.im = -sqrt(sqrt(a) * (b - a)) * sqrt((a - c) * (a - d)) / a;
    } else {
      dgamma.im *= dsdy;
    }

    break;

  case 'b':
    dsdy = 4. * y * (sqrt(b) - y * y);

    if (y >= pow(b, 1. / 4.) - 1.e-3) {
      dgamma.im = -sqrt(sqrt(b) * (b - a)) * sqrt((b - c) * (b - d)) / b;
    } else {
      dgamma.im *= dsdy;
    }

    break;

  default:
    printf("FATAL ERROR: In function dgamma_III, unknown region '%c'!", region);
    exit(1);
    break;
  }

  switch (kappa_pm) {
  case '+':
    break;

  case '-':
    dgamma.im *= -1.;
    break;

  default:
    printf("FATAL ERROR: In function dgamma_III, kappa_pm has to be '+' or "
           "'-', '%c' is unknown!",
           kappa_pm);
    exit(1);
    break;
  }

  dgamma.re *= dsdy;

  return dgamma;
}

void complex_spline_alloc(complex_spline *A, int N, int *N_spline) {
  int i;
  for (i = 0; i < N; i++) {
    A[i].re = gsl_spline_alloc(gsl_interp_cspline, N_spline[i]);
    A[i].im = gsl_spline_alloc(gsl_interp_cspline, N_spline[i]);
  }
}

void complex_spline_free(complex_spline *A, int N) {
  int i;
  for (i = 0; i < N; i++) {
    gsl_spline_free(A[i].re);
    gsl_spline_free(A[i].im);
  }
}

void cubic_spline_f01_df01(double x0, double x1, double f0, double f1,
                           double df0, double df1, double *c) {
  c[0] = -2.0 * ((f0 - f1) + 0.5 * (df0 - df1) * (x0 - x1) - df0 * (x0 - x1)) /
         ((x0 * x0 * x0 - x1 * x1 * x1) - 3.0 * x0 * x1 * (x0 - x1));
  c[1] = 0.5 * (df0 - df1) / (x0 - x1) - 1.5 * c[0] * (x0 + x1);
  c[2] = df0 - 3.0 * c[0] * x0 * x0 - 2.0 * c[1] * x0;
  c[3] = f0 - c[0] * x0 * x0 * x0 - c[1] * x0 * x0 - c[2] * x0;
}

complex R_sqrt(double s, double s0, double a, double L2) {
  double s_m_a;
  complex result;

  if (s < a) {
    s_m_a = sqrt(a - s);
    result.re = (2.0 * atanh(s_m_a / sqrt(a - s0)) - TWO_ACOTH_SQRT2) / s_m_a;
    result.im = M_PI / s_m_a;
  } else {
    s_m_a = sqrt(s - a);
    result.re = (TWO_ACOTH_SQRT2 - 2.0 * atanh(s_m_a / sqrt(L2 - a))) / s_m_a;
    result.im = M_PI / s_m_a;
  }

  return result;
}

complex R_sqrt_cubed(double s, double s0, double a, double L2, complex R12) {
  double s_m_a;
  complex result;

  if (s < a) {
    s_m_a = a - s;
    result.re =
        (2.0 * (SQRT2 / sqrt(s_m_a) - 1.0 / sqrt(a - s0)) + R12.re) / s_m_a;
    result.im = R12.im / s_m_a;

  }

  else {
    s_m_a = s - a;
    result.re =
        (2.0 * (1.0 / sqrt(L2 - a) - SQRT2 / sqrt(s_m_a)) + R12.re) / s_m_a;
    result.im = R12.im / s_m_a;
  }

  return result;
}

complex Q_sqrt(double s, double s0, double a, double L2) {
  double s_m_a;
  complex result;

  if (s < a) {
    s_m_a = sqrt(a - s);

    if (s < s0) {
      result.re = 2.0 * atanh(sqrt(a - s0) / s_m_a) / s_m_a;
    } else {
      result.re = TWO_ACOTH_SQRT2 / s_m_a;
    }

    result.im = 2.0 * atan(sqrt(L2 - a) / s_m_a) / s_m_a; //-
  }

  else {
    s_m_a = sqrt(s - a);
    result.re = -TWO_ACOTH_SQRT2 / s_m_a;
    result.im = -2.0 * atan(sqrt(a - s0) / s_m_a) / s_m_a;
  }

  return result;
}

complex Q_sqrt_cubed(double s, double s0, double a, double L2, complex Q12) {
  double s_m_a;
  complex result;

  if (s < a) {
    s_m_a = a - s;

    if (s < s0) {
      result.re = (-2.0 / sqrt(a - s0) + Q12.re) / s_m_a;
    } else {
      result.re = (-2.0 * SQRT2 / sqrt(s_m_a) + Q12.re) / s_m_a;
    }

    result.im = -(2.0 / sqrt(L2 - a) + Q12.im) / s_m_a; //(-2.0 and +Q12
  }

  else {
    s_m_a = s - a;
    result.re = (2.0 * SQRT2 / sqrt(s_m_a) + Q12.re) / s_m_a;
    result.im = (2.0 / sqrt(a - s0) - Q12.im) / s_m_a; //+(-2.0 and +Q12
  }

  return result;
}

// complex Q_sqrt_complex(complex s, double s0, double L2) {
//    double temp;
//    complex x,y,sqrt_c,log_c,result;
//
//    x.re = METAP_M_MPION_SQUARED-s.re;
//    x.im = -s.im;
//
//    x = c_sqrt(x);
//
//    sqrt_c = c_inverse(x);
//
//    temp = sqrt(METAP_M_MPION_SQUARED-s0);
//
//    x.re = sqrt_c.re*temp;
//    x.im = sqrt_c.im*temp;
//
//    x = c_two_atanh(x);
//
//    temp = sqrt(L2-METAP_M_MPION_SQUARED);
//
//    y.re = sqrt_c.re*temp;
//    y.im = sqrt_c.im*temp;
//
//    y = c_two_i_atan(y);
//
//    log_c.re = x.re-y.re;
//    log_c.im = x.im-y.im;
//
//    result.re = sqrt_c.re*log_c.re-sqrt_c.im*log_c.im;
//    result.im = sqrt_c.re*log_c.im+sqrt_c.im*log_c.re;
//
//    return result;
//}

complex Q_sqrt_complex_f(complex s, double s0, double a) {
  double temp;
  complex x, sqrt_c, result;

  x.re = a - s.re;
  x.im = -s.im;

  x = c_sqrt(x);

  sqrt_c = c_inverse(x);

  temp = sqrt(a - s0);

  x.re = sqrt_c.re * temp;
  x.im = sqrt_c.im * temp;

  x = c_two_atanh(x);

  result.re = sqrt_c.re * x.re - sqrt_c.im * x.im;
  result.im = sqrt_c.re * x.im + sqrt_c.im * x.re;

  return result;
}

complex Q_sqrt_complex_g(complex s, double a, double L2) {
  double temp;
  complex x, y, sqrt_c, result;

  x.re = a - s.re;
  x.im = -s.im;

  x = c_sqrt(x);

  sqrt_c = c_inverse(x);

  temp = sqrt(L2 - a);

  x.re = sqrt_c.re * temp;
  x.im = sqrt_c.im * temp;

  y = c_two_i_atan(x);

  x.re = y.im;
  x.im = -y.re;

  result.re = sqrt_c.re * x.re - sqrt_c.im * x.im;
  result.im = sqrt_c.re * x.im + sqrt_c.im * x.re;

  return result;
}

complex Q_sqrt_cubed_complex(complex s, double s0, double a, double L2,
                             complex Q12) {
  complex x, y, result;

  x.re = a - s.re;
  x.im = -s.im;

  x = c_inverse(x);

  y.re = -2.0 / sqrt(a - s0) + Q12.re;
  y.im = -2.0 / sqrt(L2 - a) + Q12.im;

  result.re = x.re * y.re - x.im * y.im;
  result.im = x.re * y.im + x.im * y.re;

  return result;
}

complex Q_sqrt_cubed_complex_f(complex s, double s0, double a, complex Q12_f) {
  complex x, y, result;

  x.re = a - s.re;
  x.im = -s.im;

  x = c_inverse(x);

  y.re = -2.0 / sqrt(a - s0) + Q12_f.re;
  y.im = Q12_f.im;

  result.re = x.re * y.re - x.im * y.im;
  result.im = x.re * y.im + x.im * y.re;

  return result;
}

complex Q_sqrt_cubed_complex_g(complex s, double a, double L2, complex Q12_g) {
  complex x, y, result;

  x.re = a - s.re;
  x.im = -s.im;

  x = c_inverse(x);

  y.re = -2.0 / sqrt(L2 - a) - Q12_g.re;
  y.im = -Q12_g.im;

  result.re = x.re * y.re - x.im * y.im;
  result.im = x.re * y.im + x.im * y.re;

  return result;
}

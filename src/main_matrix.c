#include "Amplitude.h"
#include "AngularAverages.h"
#include "Basic.h"
#include "Grid.h"
#include "InhomogenitySqrt.h"
#include "InhomogenitySqrtCubed.h"
#include "InputOutput.h"
#include "Omnes.h"

int main(int argn, char **argc) {
  int i, j, N[10], N_tilde, N_etapi, n, N_sin, n0, n1, n2, n0_sub, n1_sub,
      n2_sub;
  int N_M_int[4], N_M_avg[4], N_M_hat[2];
  double t, si, sf, ds, eps, eps_b, s_step, s, s0, L2, L2_tilde, L2_phase,
      s_min, s_max, k, s_step_low_high;
  double *s_cv_plus, *s_cv_minus, *s_below_cut, *phi,
      *s_tilde;       // interpolation grid
  complex *s_complex; // interpolation grid
  complex F0[4], G0[4], F1[5], G1[5], F2[4], G2[4], z, Q12, Q32, R12, R32,
      temp0, temp1, temp2;
  char sub_const[10], phase_type[100], filename[1000];
  FILE *fin, *fout;

  gsl_set_error_handler_off();

  t = clock();

  printf("Mass of decaying meson: %.3f MeV\n\n", METAP);

  // reading input file
  printf("Reading input file...\n");
  Input(argc[1], &s_step, sub_const, &n0_sub, &n1_sub, &n2_sub, phase_type);
  // printf("%.3f %s %d %d %d\n",s_step,sub_const,n0_sub,n1_sub,n2_sub);
  // printf("%s\n",filename);
  printf("Done.\n\n");

  subtraction_constant(sub_const, &n0, &n1, &n2);
  // printf("%d %d %d\n",n0,n1,n2);
  printf("Dertermine the system for subtraction constant '%s'\n\n", sub_const);

  gsl_spline *delta0;
  gsl_spline *delta1;
  gsl_spline *delta2;

  // import phases
  printf("Import scattering phase shifts...\n");
  sprintf(filename, "Input/Phases/%s_phases", phase_type);
  Phases(filename, &delta0, &delta1, &delta2, &s0, &L2_phase);
  printf("Done.\n\n");

  printf("s_th = %.6f MPION**2\n", s0);
  printf("L2_phase = %.6f MPION**2\n\n", L2_phase);

  N_tilde = 200000;
  s_tilde = (double *)malloc(N_tilde * sizeof(double));
  complex *M0_tilde = (complex *)malloc(N_tilde * sizeof(complex));
  complex *M1_tilde = (complex *)malloc(N_tilde * sizeof(complex));
  complex *M2_tilde = (complex *)malloc(N_tilde * sizeof(complex));

  eps = 1.0e-06;
  eps_b = 0.5; // treatment of dividing by 0 at b=METAP_P_MPION_SQUARED

  // import M_tilde
  printf("Import M_tilde functions...\n");
  sprintf(filename, "Input/Matrix/%s/M_tilde_%s", phase_type, sub_const);
  M_tilde_input(filename, M0_tilde, M1_tilde, M2_tilde, s_tilde, &N_M_hat[0],
                &N_M_hat[1], eps_b);
  printf("Done.\n\n");

  // import F_etapi
  printf("Import F_etapi functions...\n");
  complex_spline F0_etapi;
  complex_spline F2_etapi;
  if (str_eq_str(sub_const, "etapi") == 1) {
    double *s_F = (double *)malloc(N_tilde * sizeof(double));
    double *F0_re = (double *)malloc(N_tilde * sizeof(double));
    double *F0_im = (double *)malloc(N_tilde * sizeof(double));
    double *F2_re = (double *)malloc(N_tilde * sizeof(double));
    double *F2_im = (double *)malloc(N_tilde * sizeof(double));

    sprintf(filename, "Input/etapi_disc/%s/F_etapi", phase_type);
    F_etapi_input(filename, F0_re, F0_im, F2_re, F2_im, s_F, &N_etapi);

    F0_etapi.re = gsl_spline_alloc(gsl_interp_cspline, N_etapi);
    F0_etapi.im = gsl_spline_alloc(gsl_interp_cspline, N_etapi);
    F2_etapi.re = gsl_spline_alloc(gsl_interp_cspline, N_etapi);
    F2_etapi.im = gsl_spline_alloc(gsl_interp_cspline, N_etapi);

    gsl_spline_init(F0_etapi.re, s_F, F0_re, N_etapi);
    gsl_spline_init(F0_etapi.im, s_F, F0_im, N_etapi);
    gsl_spline_init(F2_etapi.re, s_F, F2_re, N_etapi);
    gsl_spline_init(F2_etapi.im, s_F, F2_im, N_etapi);

    free(s_F);
    free(F0_re);
    free(F0_im);
    free(F2_re);
    free(F2_im);
  }
  printf("Done.\n\n");

  // number of grid points for tilde and hat functions
  N_tilde = N_M_hat[0] + N_M_hat[1];

  L2_tilde = s_tilde[N_tilde - 1];
  printf("s_th = %.6f MPION**2\n", s0);
  printf("L2_tilde = %.6f MPION**2\n\n", L2_tilde);

  L2 = L2_tilde;

  s_max = L2;
  s_min = integration_s_plus(s_max);

  // for singularity in M_hat at a=(METAP-MPION)^2
  N_sin = (int)floor(log2(s_step / eps));
  // N_sin = 0;
  printf("N_singularity=%d %.2e\n\n", N_sin, s_step * pow(0.5, N_sin));

  s_step_low_high = 1.; // s_step for low and high energy ends
  // calculating number of grid points
  grid_point_size(N, N_sin, s_step, s_step_low_high, s_min, s_max, s0, eps_b);

  // number of grid points for angular averages
  N_M_avg[0] = N[2];
  N_M_avg[1] = N[3] + N[9];
  N_M_avg[2] = N[4];
  N_M_avg[3] = N[5] + N[6];

  // number of grid points for Omnes-functions, amplitudes and inhomogeneities
  N_M_int[0] = N_M_avg[0] + N_M_avg[1] + N_M_avg[2] + N_M_avg[3];
  N_M_int[1] = N[7];
  N_M_int[2] = N[8];
  N_M_int[3] = N[0] + N[1];

  printf("Number of grid points for the Omnes functions and amplitudes:\n");
  for (i = 0; i < 4; i++) {
    printf("N_M_int(%d)=%d\n", i + 1, N_M_int[i]);
  }
  printf("\n");

  printf("Number of grid points for angular averages:\n");
  for (i = 0; i < 4; i++) {
    printf("N_M_avg(%d)=%d\n", i + 1, N_M_avg[i]);
  }
  printf("\n");

  printf("Number of grid points for inhomogenities:\n");
  for (i = 0; i < 2; i++) {
    printf("N_M_hat(%d)=%d\n", i + 1, N_M_hat[i]);
  }
  printf("\n");

  s_below_cut = (double *)malloc(N_M_int[3] * sizeof(double));
  s_cv_plus = (double *)malloc(N_M_int[0] * sizeof(double));
  s_cv_minus = (double *)malloc(N_M_int[1] * sizeof(double));
  phi = (double *)malloc(N_M_int[2] * sizeof(double));
  s_complex = (complex *)malloc(N_M_int[2] * sizeof(complex));

  build_grid(s_cv_plus, s_cv_minus, s_below_cut, s_complex, phi, s_min, s_max,
             s0, eps, eps_b, N, N_sin);

  complex *R12_cv_plus = (complex *)malloc(N_M_int[0] * sizeof(complex));
  complex *R32_cv_plus = (complex *)malloc(N_M_int[0] * sizeof(complex));
  complex *Q12_cv_plus = (complex *)malloc(N_M_int[0] * sizeof(complex));
  complex *Q32_cv_plus = (complex *)malloc(N_M_int[0] * sizeof(complex));

  for (i = 0; i < N_M_int[0]; i++) {
    R12_cv_plus[i] = R_sqrt(s_cv_plus[i], s0, L2);
    R32_cv_plus[i] = R_sqrt_cubed(s_cv_plus[i], s0, L2, R12_cv_plus[i]);

    Q12_cv_plus[i] = Q_sqrt(s_cv_plus[i], s0, L2);
    Q32_cv_plus[i] = Q_sqrt_cubed(s_cv_plus[i], s0, L2, Q12_cv_plus[i]);
  }

  complex *R12_cv_minus = (complex *)malloc(N_M_int[1] * sizeof(complex));
  complex *R32_cv_minus = (complex *)malloc(N_M_int[1] * sizeof(complex));
  complex *Q12_cv_minus = (complex *)malloc(N_M_int[1] * sizeof(complex));
  complex *Q32_cv_minus = (complex *)malloc(N_M_int[1] * sizeof(complex));

  for (i = 0; i < N_M_int[1]; i++) {
    R12_cv_minus[i] = R_sqrt(s_cv_minus[i], s0, L2);
    R32_cv_minus[i] = R_sqrt_cubed(s_cv_minus[i], s0, L2, R12_cv_minus[i]);

    Q12_cv_minus[i] = Q_sqrt(s_cv_minus[i], s0, L2);
    Q32_cv_minus[i] = Q_sqrt_cubed(s_cv_minus[i], s0, L2, Q12_cv_minus[i]);
  }

  complex *Q12_complex_f = (complex *)malloc(N_M_int[2] * sizeof(complex));
  complex *Q32_complex_f = (complex *)malloc(N_M_int[2] * sizeof(complex));
  complex *Q12_complex_g = (complex *)malloc(N_M_int[2] * sizeof(complex));
  complex *Q32_complex_g = (complex *)malloc(N_M_int[2] * sizeof(complex));

  for (i = 0; i < N_M_int[2]; i++) {
    Q12_complex_f[i] = Q_sqrt_complex_f(s_complex[i], s0);
    Q32_complex_f[i] =
        Q_sqrt_cubed_complex_f(s_complex[i], s0, Q12_complex_f[i]);
    Q12_complex_g[i] = Q_sqrt_complex_g(s_complex[i], L2);
    Q32_complex_g[i] =
        Q_sqrt_cubed_complex_g(s_complex[i], L2, Q12_complex_g[i]);
  }

  complex *Q12_below_cut = (complex *)malloc(N_M_int[3] * sizeof(complex));
  complex *Q32_below_cut = (complex *)malloc(N_M_int[3] * sizeof(complex));

  for (i = 0; i < N_M_int[3]; i++) {
    Q12_below_cut[i] = Q_sqrt(s_below_cut[i], s0, L2);
    Q32_below_cut[i] = Q_sqrt_cubed(s_below_cut[i], s0, L2, Q12_below_cut[i]);
  }

  double *sin_delta0 = (double *)malloc(N_tilde * sizeof(double));
  double *sin_delta1 = (double *)malloc(N_tilde * sizeof(double));
  double *sin_delta2 = (double *)malloc(N_tilde * sizeof(double));

  double *abs_omnes0_dump = (double *)malloc(N_M_int[0] * sizeof(double));
  double *abs_omnes1_dump = (double *)malloc(N_M_int[0] * sizeof(double));
  double *abs_omnes2_dump = (double *)malloc(N_M_int[0] * sizeof(double));

  double *abs_omnes0 = (double *)malloc(N_tilde * sizeof(double));
  double *abs_omnes1 = (double *)malloc(N_tilde * sizeof(double));
  double *abs_omnes2 = (double *)malloc(N_tilde * sizeof(double));

  phase_shifts(delta0, delta1, delta2, sin_delta0, sin_delta1, sin_delta2,
               s_tilde, N_tilde);

  complex **omnes0 = (complex **)malloc(4 * sizeof(complex *));
  for (i = 0; i < 4; i++) {
    omnes0[i] = (complex *)malloc(N_M_int[i] * sizeof(complex));
  }

  complex **omnes1 = (complex **)malloc(4 * sizeof(complex *));
  for (i = 0; i < 4; i++) {
    omnes1[i] = (complex *)malloc(N_M_int[i] * sizeof(complex));
  }

  complex **omnes2 = (complex **)malloc(4 * sizeof(complex *));
  for (i = 0; i < 4; i++) {
    omnes2[i] = (complex *)malloc(N_M_int[i] * sizeof(complex));
  }

  complex **M0_inhom = (complex **)malloc(4 * sizeof(complex *));
  for (i = 0; i < 4; i++) {
    M0_inhom[i] = (complex *)malloc(N_M_int[i] * sizeof(complex));
  }

  complex **M1_inhom = (complex **)malloc(4 * sizeof(complex *));
  for (i = 0; i < 4; i++) {
    M1_inhom[i] = (complex *)malloc(N_M_int[i] * sizeof(complex));
  }

  complex **M2_inhom = (complex **)malloc(4 * sizeof(complex *));
  for (i = 0; i < 4; i++) {
    M2_inhom[i] = (complex *)malloc(N_M_int[i] * sizeof(complex));
  }

  n = 2;
  complex **M0_avg = (complex **)malloc(n * sizeof(complex *));
  for (i = 0; i < n; i++) {
    M0_avg[i] = (complex *)malloc(N_M_int[0] * sizeof(complex));
  }

  n = 3;
  complex **M1_avg = (complex **)malloc(n * sizeof(complex *));
  for (i = 0; i < n; i++) {
    M1_avg[i] = (complex *)malloc(N_M_int[0] * sizeof(complex));
  }

  n = 2;
  complex **M2_avg = (complex **)malloc(n * sizeof(complex *));
  for (i = 0; i < n; i++) {
    M2_avg[i] = (complex *)malloc(N_M_int[0] * sizeof(complex));
  }

  complex_spline *M0_hat = (complex_spline *)malloc(2 * sizeof(complex_spline));
  complex_spline *M1_hat = (complex_spline *)malloc(2 * sizeof(complex_spline));
  complex_spline *M2_hat = (complex_spline *)malloc(2 * sizeof(complex_spline));

  complex_spline_alloc(M0_hat, 2, N_M_hat);
  complex_spline_alloc(M1_hat, 2, N_M_hat);
  complex_spline_alloc(M2_hat, 2, N_M_hat);

  complex_spline *M0 = (complex_spline *)malloc(4 * sizeof(complex_spline));
  complex_spline *M1 = (complex_spline *)malloc(4 * sizeof(complex_spline));
  complex_spline *M2 = (complex_spline *)malloc(4 * sizeof(complex_spline));

  complex_spline_alloc(M0, 4, N_M_int);
  complex_spline_alloc(M1, 4, N_M_int);
  complex_spline_alloc(M2, 4, N_M_int);

  printf("Computing Omnes functions...\n");
  build_omnes(omnes0, abs_omnes0_dump, delta0, s0, L2_phase, 2. * M_PI,
              s_cv_plus, s_cv_minus, s_below_cut, s_complex, N_M_int);
  build_omnes(omnes1, abs_omnes1_dump, delta1, s0, L2_phase, M_PI, s_cv_plus,
              s_cv_minus, s_below_cut, s_complex, N_M_int);
  build_omnes(omnes2, abs_omnes2_dump, delta2, s0, L2_phase, 0., s_cv_plus,
              s_cv_minus, s_below_cut, s_complex, N_M_int);

  abs_omnes_cv_plus(abs_omnes0, delta0, s0, L2_phase, 2. * M_PI, s_tilde,
                    N_tilde);
  abs_omnes_cv_plus(abs_omnes1, delta1, s0, L2_phase, M_PI, s_tilde, N_tilde);
  abs_omnes_cv_plus(abs_omnes2, delta2, s0, L2_phase, 0., s_tilde, N_tilde);
  printf("Done.\n\n");

  printf("Building inhomogenity integrands...\n");
  build_inhomogenity_integrand_Matrix(
      M0_hat, M1_hat, M2_hat, F0, G0, F1, G1, F2, G2, M0_tilde, M1_tilde,
      M2_tilde, sin_delta0, sin_delta1, sin_delta2, abs_omnes0, abs_omnes1,
      abs_omnes2, s_tilde, s0, L2, N_M_hat, n0_sub, n1_sub, n2_sub);
  printf("Done.\n\n");

  printf("Computing inhomogenities...\n");
  build_inhomogenity_sqrt(
      M0_inhom, M0_hat, F0, G0, s0, L2, s_cv_plus, s_cv_minus, s_below_cut,
      s_complex, R12_cv_plus, Q12_cv_plus, R12_cv_minus, Q12_cv_minus,
      Q12_below_cut, Q12_complex_f, Q12_complex_g, N_M_int, n0_sub);
  build_inhomogenity_sqrt_cubed(
      M1_inhom, M1_hat, F1, G1, s0, L2, s_cv_plus, s_cv_minus, s_below_cut,
      s_complex, R32_cv_plus, Q12_cv_plus, Q32_cv_plus, R32_cv_minus,
      Q12_cv_minus, Q32_cv_minus, Q12_below_cut, Q32_below_cut, Q12_complex_f,
      Q32_complex_f, Q12_complex_g, Q32_complex_g, N_M_int, n1_sub);
  build_inhomogenity_sqrt(
      M2_inhom, M2_hat, F2, G2, s0, L2, s_cv_plus, s_cv_minus, s_below_cut,
      s_complex, R12_cv_plus, Q12_cv_plus, R12_cv_minus, Q12_cv_minus,
      Q12_below_cut, Q12_complex_f, Q12_complex_g, N_M_int, n2_sub);
  printf("Done.\n\n");

  gsl_interp_accel *acc = gsl_interp_accel_alloc();
  if (str_eq_str(sub_const, "etapi") == 1) {
    for (i = 0; i < N_M_int[0]; i++) {
      M0_inhom[0][i].re += gsl_spline_eval(F0_etapi.re, s_cv_plus[i], acc);
      M0_inhom[0][i].im += gsl_spline_eval(F0_etapi.im, s_cv_plus[i], acc);
      M2_inhom[0][i].re += gsl_spline_eval(F2_etapi.re, s_cv_plus[i], acc);
      M2_inhom[0][i].im += gsl_spline_eval(F2_etapi.im, s_cv_plus[i], acc);
    }
    for (i = 0; i < N_M_int[3]; i++) {
      M0_inhom[3][i].re += gsl_spline_eval(F0_etapi.re, s_below_cut[i], acc);
      M0_inhom[3][i].im += gsl_spline_eval(F0_etapi.im, s_below_cut[i], acc);
      M2_inhom[3][i].re += gsl_spline_eval(F2_etapi.re, s_below_cut[i], acc);
      M2_inhom[3][i].im += gsl_spline_eval(F2_etapi.im, s_below_cut[i], acc);
    }
  }
  gsl_interp_accel_free(acc);

  printf("Computing Amplitudes...\n");
  build_amplitude(M0, omnes0, M0_inhom, s_cv_plus, s_cv_minus, s_below_cut,
                  s_complex, phi, s0, L2, s_min, N_M_int, n0);
  build_amplitude(M1, omnes1, M1_inhom, s_cv_plus, s_cv_minus, s_below_cut,
                  s_complex, phi, s0, L2, s_min, N_M_int, n1);
  build_amplitude(M2, omnes2, M2_inhom, s_cv_plus, s_cv_minus, s_below_cut,
                  s_complex, phi, s0, L2, s_min, N_M_int, n2);
  printf("Done.\n\n");

  printf("Computing angular averages...\n");
  for (n = 0; n < 2; n++) {
    build_M_avg(M0, M0_avg[n], s_cv_plus, N_M_avg, n);
    build_M_avg(M1, M1_avg[n], s_cv_plus, N_M_avg, n);
    build_M_avg(M2, M2_avg[n], s_cv_plus, N_M_avg, n);
  }
  build_M_avg(M1, M1_avg[n], s_cv_plus, N_M_avg, n);
  printf("Done.\n\n");

  printf("Writing output...\n");

  Omnes_output("OmnesFunctions", omnes0, omnes1, omnes2, s_below_cut, s_cv_plus,
               N_M_int[3], N_M_int[0]);

  sprintf(filename, "M_tilde_%s", sub_const);
  Inhomogenities_output(filename, M0_avg, M1_avg, M2_avg, s_cv_plus,
                        N_M_int[0]);

  sprintf(filename, "M_int_%s", sub_const);
  Dispersion_integral_output(filename, M0_inhom, M1_inhom, M2_inhom,
                             s_below_cut, s_cv_plus, N_M_int[3], N_M_int[0]);

  sprintf(filename, "M_%s", sub_const);
  Amplitudes_output(filename, M0, M1, M2, s_below_cut, s_cv_plus, N_M_int[3],
                    N_M_int[0]);

  printf("Done.\n\n");

  free(s_cv_plus);
  free(s_cv_minus);
  free(s_complex);
  free(phi);
  free(s_below_cut);

  free(R12_cv_plus);
  free(R32_cv_plus);
  free(Q12_cv_plus);
  free(Q32_cv_plus);

  free(R12_cv_minus);
  free(R32_cv_minus);
  free(Q12_cv_minus);
  free(Q32_cv_minus);

  free(Q12_complex_f);
  free(Q32_complex_f);
  free(Q12_complex_g);
  free(Q32_complex_g);

  free(Q12_below_cut);
  free(Q32_below_cut);

  free(sin_delta0);
  free(sin_delta1);
  free(sin_delta2);
  free(abs_omnes0);
  free(abs_omnes1);
  free(abs_omnes2);
  free(abs_omnes0_dump);
  free(abs_omnes1_dump);
  free(abs_omnes2_dump);

  for (i = 0; i < 4; i++) {
    free(omnes0[i]);
  }
  free(omnes0);

  for (i = 0; i < 4; i++) {
    free(omnes1[i]);
  }
  free(omnes1);

  for (i = 0; i < 4; i++) {
    free(omnes2[i]);
  }
  free(omnes2);

  for (i = 0; i < 4; i++) {
    free(M0_inhom[i]);
  }
  free(M0_inhom);

  for (i = 0; i < 4; i++) {
    free(M1_inhom[i]);
  }
  free(M1_inhom);

  for (i = 0; i < 4; i++) {
    free(M2_inhom[i]);
  }
  free(M2_inhom);

  n = 2;
  for (i = 0; i < n; i++) {
    free(M0_avg[i]);
  }
  free(M0_avg);

  n = 3;
  for (i = 0; i < n; i++) {
    free(M1_avg[i]);
  }
  free(M1_avg);

  n = 2;
  for (i = 0; i < n; i++) {
    free(M2_avg[i]);
  }
  free(M2_avg);

  complex_spline_free(M0_hat, 2);
  complex_spline_free(M1_hat, 2);
  complex_spline_free(M2_hat, 2);
  free(M0_hat);
  free(M1_hat);
  free(M2_hat);

  complex_spline_free(M0, 4);
  complex_spline_free(M1, 4);
  complex_spline_free(M2, 4);
  free(M0);
  free(M1);
  free(M2);

  gsl_spline_free(delta0);
  gsl_spline_free(delta1);
  gsl_spline_free(delta2);

  t = (clock() - t) / CLOCKS_PER_SEC;

  printf("\ncalculation time = %.3fs\n\n", t);

  return 0;
}

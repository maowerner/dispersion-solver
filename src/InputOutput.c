#include "InputOutput.h"
#include "Basic.h"
#include "Singularity.h"

// reading input file
void input(char *filename, double *s_step, double *L2, int *n0_sub, int *n1_sub,
           int *n2_sub) {
  int i;
  char charged_or_neutral;
  FILE *fin;

  if ((fin = fopen(filename, "r")) == NULL) {
    printf("FATAL ERROR: Bad input, no input file.\n");
    fclose(fin);
    exit(1);
  }

  fscanf(fin, "%lf %lf %d %c", s_step, L2, &N_ITER, &charged_or_neutral);
  fscanf(fin, "%d %d %d", n0_sub, n1_sub, n2_sub);
  fscanf(fin, "%d", &N_SUB_CONST);

  switch (charged_or_neutral) {
  case 'C':
    pi0pi0 = false;
    printf("Calculate the system for charged pions.\n");
    break;

  case 'N':
    pi0pi0 = true;
    printf("Calculate the system for neutral pions.\n");
    break;

  default:
    printf("\nFATAL ERROR: In input file. Unknown pipi System, choose 'C' for "
           "pi+pi- or 'N' for pi0pi0 system.\n");
    exit(1);
    break;
  }

  SUB_CONST = (char **)malloc(N_SUB_CONST * sizeof(char *));
  for (i = 0; i < N_SUB_CONST; i++) {
    SUB_CONST[i] = (char *)malloc(2 * sizeof(char));
  }

  for (i = 0; i < N_SUB_CONST; i++) {
    fscanf(fin, "%s", SUB_CONST[i]);
  }

  fscanf(fin, "%s %lf %lf %lf", &file_delta0, &delta0_params.L2,
         &delta0_params.n, &delta0_params.m);
  fscanf(fin, "%s %lf %lf %lf", &file_delta1, &delta1_params.L2,
         &delta1_params.n, &delta1_params.m);
  fscanf(fin, "%s %lf %lf %lf", &file_delta2, &delta2_params.L2,
         &delta2_params.n, &delta2_params.m);
  fscanf(fin, "%s %lf %lf %lf", &file_delta_etapi, &delta_etapi_params.L2,
         &delta_etapi_params.n, &delta_etapi_params.m);

  fclose(fin);
}

// read in phase shift input from txt
void phase_input(char *filename, gsl_spline **delta_spline) {
  int i, N, n;
  double *s, *phase;
  FILE *fin;

  if ((fin = fopen(filename, "r")) == NULL) {
    printf("FATAL ERROR: Bad input, file '%s' does not exist.\n", filename);
    fclose(fin);
    exit(1);
  }

  fscanf(fin, "%d", &N);

  if ((s = (double *)malloc(N * sizeof(double))) == NULL) {
    printf("FATAL ERROR: In memory allocation 's'.\n");
    exit(1);
  }
  if ((phase = (double *)malloc(N * sizeof(double))) == NULL) {
    printf("FATAL ERROR: In memory allocation 'phase'.\n");
    exit(1);
  }

  for (i = 0; i < N; i++) {
    fscanf(fin, "%lf %lf", &s[i], &phase[i]);
    if (i > 0) {
      if (s[i - 1] >= s[i]) {
        printf("FATAL ERROR: while reading phase input of file '%s' "
               "s[i-1]>=s[i] at i=%d!\n\n",
               filename, i);
        exit(1);
      }
    }
  }

  fclose(fin);

  *delta_spline = gsl_spline_alloc(gsl_interp_cspline, N);

  gsl_spline_init(*delta_spline, s, phase, N);

  free(s);
  free(phase);
}

void M_tilde_input(char *filename, complex *M0_tilde, complex *M1_tilde,
                   complex *M2_tilde, double *s, int *N_low, int *N_high,
                   double eps_b) {
  int i, N, n;
  double si;
  FILE *fin;

  if ((fin = fopen(filename, "r")) == NULL) {
    printf("FATAL ERROR: Bad input, file '%s' does not exist.\n", filename);
    fclose(fin);
    exit(1);
  }

  fscanf(fin, "%d", &N);
  N--; // ignoring s[0]=s0

  //    if ((*s = (double *)malloc(N*sizeof(double)))==NULL) {
  //        printf("FATAL ERROR: In memory allocation 's'.\n");
  //        exit(1);
  //    }
  //    if ((*M0_tilde = (complex *)malloc(N*sizeof(complex)))==NULL) {
  //        printf("FATAL ERROR: In memory allocation 'M0_tilde'.\n");
  //        exit(1);
  //    }
  //    if ((*M1_tilde = (complex *)malloc(N*sizeof(complex)))==NULL) {
  //        printf("FATAL ERROR: In memory allocation 'M1_re'.\n");
  //        exit(1);
  //    }
  //    if ((*M2_tilde = (complex *)malloc(N*sizeof(complex)))==NULL) {
  //        printf("FATAL ERROR: In memory allocation 'M2_re'.\n");
  //        exit(1);
  //    }

  printf("%d\n", N);
  fscanf(fin, "%lf%*lf%*lf%*lf%*lf%*lf%*lf", &si); // ignoring s[0]=s0

  (*N_low) = 0;
  (*N_high) = 0;
  for (i = 0; si < METAP_P_MPION_SQUARED - eps_b; i++) {
    fscanf(fin, "%lf%lf%lf%lf%lf%lf%lf", &si, &M0_tilde[i].re, &M0_tilde[i].im,
           &M1_tilde[i].re, &M1_tilde[i].im, &M2_tilde[i].re, &M2_tilde[i].im);
    s[i] = si;
    if (s[i] < METAP_M_MPION_SQUARED) {
      (*N_low)++;
    }
    if (s[i] > METAP_M_MPION_SQUARED) {
      (*N_high)++;
    }
  }
  n = 0;
  for (; si < METAP_P_MPION_SQUARED + eps_b; i++) {
    fscanf(fin, "%lf%*lf%*lf%*lf%*lf%*lf%*lf", &si);
    n++;
  }

  for (; i < N; i++) {
    fscanf(fin, "%lf%lf%lf%lf%lf%lf%lf", &s[i - n], &M0_tilde[i - n].re,
           &M0_tilde[i - n].im, &M1_tilde[i - n].re, &M1_tilde[i - n].im,
           &M2_tilde[i - n].re, &M2_tilde[i - n].im);
    (*N_high)++;
  }

  fclose(fin);
}

void F_etapi_input(char *filename, double *F0_re, double *F0_im, double *F2_re,
                   double *F2_im, double *s, int *N) {
  int i;
  FILE *fin;

  if ((fin = fopen(filename, "r")) == NULL) {
    printf("FATAL ERROR: Bad input, file '%s' does not exist.\n", filename);
    fclose(fin);
    exit(1);
  }

  fscanf(fin, "%d", &N[0]);
  for (i = 0; i < N[0]; i++) {
    fscanf(fin, "%lf%lf%lf%*lf%*lf%lf%lf", &s[i], &F0_re[i], &F0_im[i],
           &F2_re[i], &F2_im[i]);
  }
  fclose(fin);
}

void Partial_wave_eta_3pi_input(char *filename, complex_spline f, double L2,
                                int n) {
  int i;
  double *s, *f_re, *f_im;
  FILE *in;

  s = (double *)malloc(n * sizeof(double));
  f_re = (double *)malloc(n * sizeof(double));
  f_im = (double *)malloc(n * sizeof(double));

  if ((in = fopen(filename, "r")) == NULL) {
    printf("FATAL ERROR: Datafile '%s' does not exist!\n\n", filename);
    exit(1);
  }
  for (i = 0; i < n; i++) {
    fscanf(in, "%lf%lf%lf", &s[i], &f_re[i], &f_im[i]);
  }

  s[0] = s_etapi;
  s[n - 1] = L2;
  gsl_spline_init(f.re, s, f_re, n);
  gsl_spline_init(f.im, s, f_im, n);

  free(s);
  free(f_re);
  free(f_im);

  fclose(in);
}

void Partial_wave_etap_eta2pi_input(char *filename, complex_spline M,
                                    complex_spline *M_tilde, double L2, int n,
                                    int m) {
  int i, j;
  double *s, *M_re, *M_im, *M_tilde_re, *M_tilde_im;
  FILE *in;

  s = (double *)malloc(n * sizeof(double));
  M_re = (double *)malloc(n * sizeof(double));
  M_im = (double *)malloc(n * sizeof(double));
  M_tilde_re = (double *)malloc(n * sizeof(double));
  M_tilde_im = (double *)malloc(n * sizeof(double));

  if ((in = fopen(filename, "r")) == NULL) {
    printf("FATAL ERROR: Datafile '%s' does not exist!\n\n", filename);
    exit(1);
  }

  for (i = 0; i < n; i++) {
    fscanf(in, "%lf%lf%lf%lf%lf", &s[i], &M_re[i], &M_im[i], &M_tilde_re[i],
           &M_tilde_im[i]);
  }
  // s[0] = s_etapi;
  s[n - 1] = L2;

  gsl_spline_init(M.re, s, M_re, n);
  gsl_spline_init(M.im, s, M_im, n);

  s[m - 1] = METAP_M_MPION_SQUARED;
  gsl_spline_init(M_tilde[0].re, s, M_tilde_re, m);
  gsl_spline_init(M_tilde[0].im, s, M_tilde_im, m);

  s[m] = METAP_M_MPION_SQUARED;
  gsl_spline_init(M_tilde[1].re, s + m, M_tilde_re + m, n - m);
  gsl_spline_init(M_tilde[1].im, s + m, M_tilde_im + m, n - m);

  free(s);
  free(M_re);
  free(M_im);
  free(M_tilde_re);
  free(M_tilde_im);

  fclose(in);
}

void Output_test(char *filename, complex_spline f0, complex_spline f2,
                 complex_spline M, complex_spline *M_tilde, double ds,
                 double L2) {
  int i;
  double s, s0, eps;

  FILE *out;

  gsl_interp_accel *acc0 = gsl_interp_accel_alloc();
  gsl_interp_accel *acc1 = gsl_interp_accel_alloc();
  gsl_interp_accel *acc2 = gsl_interp_accel_alloc();
  gsl_interp_accel *acc3 = gsl_interp_accel_alloc();

  if ((out = fopen(filename, "w")) == NULL) {
    printf("FATAL ERROR: Datafile '%s' does not exist!\n\n", filename);
    exit(1);
  }

  eps = 1.e-15;
  s0 = s_etapi;
  s = s0;
  for (i = 0; s <= L2; i++) {
    if (s < METAP_M_MPION_SQUARED) {
      fprintf(
          out,
          "%+.10e %+.10e %+.10e %+.10e %+.10e %+.10e %+.10e %+.10e %+.10e\n", s,
          gsl_spline_eval(f0.re, s, acc0), gsl_spline_eval(f0.im, s, acc0),
          gsl_spline_eval(f2.re, s, acc0), gsl_spline_eval(f2.im, s, acc0),
          gsl_spline_eval(M.re, s, acc1), gsl_spline_eval(M.im, s, acc1),
          gsl_spline_eval(M_tilde[0].re, s, acc2),
          gsl_spline_eval(M_tilde[0].im, s, acc2));
    } else if (s - ds < METAP_M_MPION_SQUARED && s > METAP_M_MPION_SQUARED) {
      s = METAP_M_MPION_SQUARED;
      fprintf(
          out,
          "%+.10e %+.10e %+.10e %+.10e %+.10e %+.10e %+.10e %+.10e %+.10e\n", s,
          gsl_spline_eval(f0.re, s, acc0), gsl_spline_eval(f0.im, s, acc0),
          gsl_spline_eval(f2.re, s, acc0), gsl_spline_eval(f2.im, s, acc0),
          gsl_spline_eval(M.re, s, acc1), gsl_spline_eval(M.im, s, acc1),
          gsl_spline_eval(M_tilde[0].re, s, acc2),
          gsl_spline_eval(M_tilde[0].im, s, acc2));
      fprintf(
          out,
          "%+.10e %+.10e %+.10e %+.10e %+.10e %+.10e %+.10e %+.10e %+.10e\n", s,
          gsl_spline_eval(f0.re, s, acc0), gsl_spline_eval(f0.im, s, acc0),
          gsl_spline_eval(f2.re, s, acc0), gsl_spline_eval(f2.im, s, acc0),
          gsl_spline_eval(M.re, s, acc1), gsl_spline_eval(M.im, s, acc1),
          gsl_spline_eval(M_tilde[1].re, s, acc3),
          gsl_spline_eval(M_tilde[1].im, s, acc3));
    } else {
      fprintf(
          out,
          "%+.10e %+.10e %+.10e %+.10e %+.10e %+.10e %+.10e %+.10e %+.10e\n", s,
          gsl_spline_eval(f0.re, s, acc0), gsl_spline_eval(f0.im, s, acc0),
          gsl_spline_eval(f2.re, s, acc0), gsl_spline_eval(f2.im, s, acc0),
          gsl_spline_eval(M.re, s, acc1), gsl_spline_eval(M.im, s, acc1),
          gsl_spline_eval(M_tilde[1].re, s, acc3),
          gsl_spline_eval(M_tilde[1].im, s, acc3));
    }

    s = s0 + ds * i;
  }

  gsl_interp_accel_free(acc0);
  gsl_interp_accel_free(acc1);
  gsl_interp_accel_free(acc2);
  gsl_interp_accel_free(acc3);

  fclose(out);
}

void Omnes_output(char *filename, complex **omnes0, complex **omnes1,
                  complex **omnes2, double *s_below_cut, double *s_cv_plus,
                  int N_below_cut, int N_cv_plus) {
  int i;
  double s;
  FILE *fout;

  fout = fopen(filename, "w");
  for (i = 0; i < N_below_cut; i++) {
    s = s_below_cut[i];
    fprintf(fout, "%+.15e %+.15e %+.15e %+.15e %+.15e %+.15e %+.15e\n", s,
            omnes0[2][i].re, omnes0[2][i].im, omnes1[2][i].re, omnes1[2][i].im,
            omnes2[2][i].re, omnes2[2][i].im);
  }
  for (i = 0; i < N_cv_plus; i++) {
    s = s_cv_plus[i];
    fprintf(fout, "%+.15e %+.15e %+.15e %+.15e %+.15e %+.15e %+.15e\n", s,
            omnes0[0][i].re, omnes0[0][i].im, omnes1[0][i].re, omnes1[0][i].im,
            omnes2[0][i].re, omnes2[0][i].im);
  }
  fclose(fout);
}

void Inhomogenities_output(char *filename, complex **M0_avg, complex **M1_avg,
                           complex **M2_avg, double *s, int N) {
  int i;
  double sigma, k;
  complex M0_hat, M1_hat, M2_hat;
  FILE *fout;

  fout = fopen(filename, "w");
  for (i = 0; i < N; i++) {
    sigma = s[i] - (METAP_SQUARED / 3. + 1.);
    k = 1.; // 1./sqrt(fabs((1.-4./s[i])*(s[i]-METAP_M_MPION_SQUARED)*(s[i]-METAP_P_MPION_SQUARED)));

    M0_hat.re = (2. / 3. * M0_avg[0][i].re + 2. * sigma * M1_avg[0][i].re +
                 4. / 3. * M1_avg[1][i].re + 20. / 9. * M2_avg[0][i].re) *
                k;
    M0_hat.im = (2. / 3. * M0_avg[0][i].im + 2. * sigma * M1_avg[0][i].im +
                 4. / 3. * M1_avg[1][i].im + 20. / 9. * M2_avg[0][i].im) *
                k;

    M1_hat.re = (6. * M0_avg[1][i].re + 9. * sigma * M1_avg[1][i].re +
                 6. * M1_avg[2][i].re - 10. * M2_avg[1][i].re) *
                k * k * k;
    M1_hat.im = (6. * M0_avg[1][i].im + 9. * sigma * M1_avg[1][i].im +
                 6. * M1_avg[2][i].im - 10. * M2_avg[1][i].im) *
                k * k * k;

    M2_hat.re = (M0_avg[0][i].re - 3. / 2. * sigma * M1_avg[0][i].re -
                 M1_avg[1][i].re + 1. / 3. * M2_avg[0][i].re) *
                k;
    M2_hat.im = (M0_avg[0][i].im - 3. / 2. * sigma * M1_avg[0][i].im -
                 M1_avg[1][i].im + 1. / 3. * M2_avg[0][i].im) *
                k;

    fprintf(fout, "%+.15e %+.15e %+.15e %+.15e %+.15e %+.15e %+.15e\n", s[i],
            M0_hat.re, M0_hat.im, M1_hat.re, M1_hat.im, M2_hat.re, M2_hat.im);

    // fprintf(fout,"%+.15e %+.15e %+.15e %+.15e %+.15e %+.15e %+.15e %+.15e
    // %+.15e %+.15e %+.15e %+.15e %+.15e %+.15e
    // %+.15e\n",s[i],M0_avg[0][i].re,M0_avg[0][i].im,M0_avg[1][i].re,M0_avg[1][i].im,M1_avg[0][i].re,M1_avg[0][i].im,M1_avg[1][i].re,M1_avg[1][i].im,M1_avg[2][i].re,M1_avg[2][i].im,M2_avg[0][i].re,M2_avg[0][i].im,M2_avg[1][i].re,M2_avg[1][i].im);

    //        if (s[i]<METAP_M_MPION_SQUARED) {
    //            fprintf(fout,"%+.15e %+.15e %+.15e %+.15e %+.15e %+.15e
    //            %+.15e\n",s[i],M0_hat.re,M0_hat.im,M1_hat.re,M1_hat.im,M2_hat.re,M2_hat.im);
    //        }
    //        else if (s[i]<METAP_P_MPION_SQUARED) {
    //            fprintf(fout,"%+.15e %+.15e %+.15e %+.15e %+.15e %+.15e
    //            %+.15e\n",s[i],M0_hat.im,-M0_hat.re,-M1_hat.im,M1_hat.re,M2_hat.im,-M2_hat.re);
    //        }
    //        else {
    //            fprintf(fout,"%+.15e %+.15e %+.15e %+.15e %+.15e %+.15e
    //            %+.15e\n",s[i],-M0_hat.re,-M0_hat.im,-M1_hat.re,-M1_hat.im,-M2_hat.re,-M2_hat.im);
    //        }
  }
  fclose(fout);
}

void Dispersion_integral_output(char *filename, complex **M0_inhom,
                                complex **M1_inhom, complex **M2_inhom,
                                double *s_below_cut, double *s_cv_plus,
                                int N_below_cut, int N_cv_plus) {
  int i;
  double s;
  FILE *fout;

  fout = fopen(filename, "w");
  for (i = 0; i < N_below_cut; i++) {
    s = s_below_cut[i];
    fprintf(fout, "%+.15e %+.15e %+.15e %+.15e %+.15e %+.15e %+.15e\n", s,
            M0_inhom[2][i].re, M0_inhom[2][i].im, M1_inhom[2][i].re,
            M1_inhom[2][i].im, M2_inhom[2][i].re, M2_inhom[2][i].im);
  }
  for (i = 0; i < N_cv_plus; i++) {
    s = s_cv_plus[i];
    fprintf(fout, "%+.15e %+.15e %+.15e %+.15e %+.15e %+.15e %+.15e\n", s,
            M0_inhom[0][i].re, M0_inhom[0][i].im, M1_inhom[0][i].re,
            M1_inhom[0][i].im, M2_inhom[0][i].re, M2_inhom[0][i].im);
  }
  fclose(fout);
}

void Amplitudes_output(char *filename, complex_spline *M0, complex_spline *M1,
                       complex_spline *M2, double *s_below_cut,
                       double *s_cv_plus, int N_below_cut, int N_cv_plus) {
  int i;
  double s;
  FILE *fout;

  gsl_interp_accel *acc0 = gsl_interp_accel_alloc();
  gsl_interp_accel *acc3 = gsl_interp_accel_alloc();

  fout = fopen(filename, "w");
  for (i = 0; i < N_below_cut; i++) {
    s = s_below_cut[i];
    fprintf(
        fout, "%+.15e %+.15e %+.15e %+.15e %+.15e %+.15e %+.15e\n", s,
        gsl_spline_eval(M0[2].re, s, acc3), gsl_spline_eval(M0[2].im, s, acc3),
        gsl_spline_eval(M1[2].re, s, acc3), gsl_spline_eval(M1[2].im, s, acc3),
        gsl_spline_eval(M2[2].re, s, acc3), gsl_spline_eval(M2[2].im, s, acc3));
  }
  for (i = 0; i < N_cv_plus; i++) {
    s = s_cv_plus[i];
    fprintf(
        fout, "%+.15e %+.15e %+.15e %+.15e %+.15e %+.15e %+.15e\n", s,
        gsl_spline_eval(M0[0].re, s, acc0), gsl_spline_eval(M0[0].im, s, acc0),
        gsl_spline_eval(M1[0].re, s, acc0), gsl_spline_eval(M1[0].im, s, acc0),
        gsl_spline_eval(M2[0].re, s, acc0), gsl_spline_eval(M2[0].im, s, acc0));
  }
  fclose(fout);

  gsl_interp_accel_free(acc0);
  gsl_interp_accel_free(acc3);
}

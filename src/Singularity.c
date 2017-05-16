#include "Singularity.h"
#include "Basic.h"

void singularity_sqrt(complex *f, double *s, complex *c) {
  double s0, s1, s2, s3;
  s0 = s[0];
  s1 = s[1];
  s2 = s[2];
  s3 = s[3];

  c[0].re =
      (f[3].re * s0 * (s0 - s1) * s1 * (s0 - s2) * (s1 - s2) * s2 +
       (-f[2].re * s0 * (s0 - s1) * s1 * (s0 - s3) * (s1 - s3) +
        s2 * (f[1].re * s0 * (s0 - s2) * (s0 - s3) -
              f[0].re * s1 * (s1 - s2) * (s1 - s3)) *
            (s2 - s3)) *
           s3) /
      ((s0 - s1) * (s0 - s2) * (s1 - s2) * (s0 - s3) * (s1 - s3) * (s2 - s3));

  c[1].re =
      (f[3].re * (s1 * s1 * s2 * s2 * (-s1 + s2) +
                  s0 * s0 * s0 * (-s1 * s1 + s2 * s2) +
                  s0 * s0 * (s1 * s1 * s1 - s2 * s2 * s2)) +
       f[2].re *
           (s1 * s1 * (s1 - s3) * s3 * s3 + s0 * s0 * s0 * (s1 * s1 - s3 * s3) +
            s0 * s0 * (-s1 * s1 * s1 + s3 * s3 * s3)) -
       (s2 - s3) *
           (f[1].re * (s0 - s2) * (s0 - s3) * (s2 * s3 + s0 * (s2 + s3)) -
            f[0].re * (s1 - s2) * (s1 - s3) * (s2 * s3 + s1 * (s2 + s3)))) /
      ((s0 - s1) * (s0 - s2) * (s1 - s2) * (s0 - s3) * (s1 - s3) * (s2 - s3));

  c[2].re =
      (f[1].re * s0 * s0 * s0 * s2 - f[0].re * s1 * s1 * s1 * s2 -
       f[1].re * s0 * s2 * s2 * s2 + f[0].re * s1 * s2 * s2 * s2 +
       f[3].re * (s0 * s0 * s0 * (s1 - s2) + s1 * s2 * (s1 * s1 - s2 * s2) +
                  s0 * (-s1 * s1 * s1 + s2 * s2 * s2)) -
       f[1].re * s0 * s0 * s0 * s3 + f[0].re * s1 * s1 * s1 * s3 -
       f[0].re * s2 * s2 * s2 * s3 + f[1].re * s2 * s2 * s2 * s3 +
       f[1].re * s0 * s3 * s3 * s3 - f[0].re * s1 * s3 * s3 * s3 +
       f[0].re * s2 * s3 * s3 * s3 - f[1].re * s2 * s3 * s3 * s3 +
       f[2].re *
           (-s1 * s1 * s1 * s3 + s1 * s3 * s3 * s3 + s0 * s0 * s0 * (-s1 + s3) +
            s0 * (s1 * s1 * s1 - s3 * s3 * s3))) /
      ((s0 - s1) * (s0 - s2) * (s1 - s2) * (s0 - s3) * (s1 - s3) * (s2 - s3));

  c[3].re =
      (-f[3].re * (s0 - s1) * (s0 - s2) * (s1 - s2) +
       f[2].re * (s0 - s1) * (s0 - s3) * (s1 - s3) -
       (f[1].re * (s0 - s2) * (s0 - s3) - f[0].re * (s1 - s2) * (s1 - s3)) *
           (s2 - s3)) /
      ((s0 - s1) * (s0 - s2) * (s1 - s2) * (s0 - s3) * (s1 - s3) * (s2 - s3));

  c[0].im =
      (f[3].im * s0 * (s0 - s1) * s1 * (s0 - s2) * (s1 - s2) * s2 +
       (-f[2].im * s0 * (s0 - s1) * s1 * (s0 - s3) * (s1 - s3) +
        s2 * (f[1].im * s0 * (s0 - s2) * (s0 - s3) -
              f[0].im * s1 * (s1 - s2) * (s1 - s3)) *
            (s2 - s3)) *
           s3) /
      ((s0 - s1) * (s0 - s2) * (s1 - s2) * (s0 - s3) * (s1 - s3) * (s2 - s3));

  c[1].im =
      (f[3].im * (s1 * s1 * s2 * s2 * (-s1 + s2) +
                  s0 * s0 * s0 * (-s1 * s1 + s2 * s2) +
                  s0 * s0 * (s1 * s1 * s1 - s2 * s2 * s2)) +
       f[2].im *
           (s1 * s1 * (s1 - s3) * s3 * s3 + s0 * s0 * s0 * (s1 * s1 - s3 * s3) +
            s0 * s0 * (-s1 * s1 * s1 + s3 * s3 * s3)) -
       (s2 - s3) *
           (f[1].im * (s0 - s2) * (s0 - s3) * (s2 * s3 + s0 * (s2 + s3)) -
            f[0].im * (s1 - s2) * (s1 - s3) * (s2 * s3 + s1 * (s2 + s3)))) /
      ((s0 - s1) * (s0 - s2) * (s1 - s2) * (s0 - s3) * (s1 - s3) * (s2 - s3));

  c[2].im =
      (f[1].im * s0 * s0 * s0 * s2 - f[0].im * s1 * s1 * s1 * s2 -
       f[1].im * s0 * s2 * s2 * s2 + f[0].im * s1 * s2 * s2 * s2 +
       f[3].im * (s0 * s0 * s0 * (s1 - s2) + s1 * s2 * (s1 * s1 - s2 * s2) +
                  s0 * (-s1 * s1 * s1 + s2 * s2 * s2)) -
       f[1].im * s0 * s0 * s0 * s3 + f[0].im * s1 * s1 * s1 * s3 -
       f[0].im * s2 * s2 * s2 * s3 + f[1].im * s2 * s2 * s2 * s3 +
       f[1].im * s0 * s3 * s3 * s3 - f[0].im * s1 * s3 * s3 * s3 +
       f[0].im * s2 * s3 * s3 * s3 - f[1].im * s2 * s3 * s3 * s3 +
       f[2].im *
           (-s1 * s1 * s1 * s3 + s1 * s3 * s3 * s3 + s0 * s0 * s0 * (-s1 + s3) +
            s0 * (s1 * s1 * s1 - s3 * s3 * s3))) /
      ((s0 - s1) * (s0 - s2) * (s1 - s2) * (s0 - s3) * (s1 - s3) * (s2 - s3));

  c[3].im =
      (-f[3].im * (s0 - s1) * (s0 - s2) * (s1 - s2) +
       f[2].im * (s0 - s1) * (s0 - s3) * (s1 - s3) -
       (f[1].im * (s0 - s2) * (s0 - s3) - f[0].im * (s1 - s2) * (s1 - s3)) *
           (s2 - s3)) /
      ((s0 - s1) * (s0 - s2) * (s1 - s2) * (s0 - s3) * (s1 - s3) * (s2 - s3));
}

void singularity_sqrt_cubed(complex *f, double *s, complex *c) {
  double s0, s1, s2, s3, s4, A, B, C, D, E;
  s0 = s[0];
  s1 = s[1];
  s2 = s[2];
  s3 = s[3];
  s4 = s[4];

  A = f[0].re;
  B = f[1].re;
  C = f[2].re;
  D = f[3].re;
  E = f[4].re;

  c[0].re = (E * s0 * s0 * (s0 - s1) * s1 * s1 * (s0 - s2) * (s1 - s2) * s2 *
                 s2 * (s0 - s3) * (s1 - s3) * (s2 - s3) * s3 * s3 +
             (-D * s0 * s0 * (s0 - s1) * s1 * s1 * (s0 - s2) * (s1 - s2) * s2 *
                  s2 * (s0 - s4) * (s1 - s4) * (s2 - s4) +
              s3 * s3 * (C * s0 * s0 * (s0 - s1) * s1 * s1 * (s0 - s3) *
                             (s1 - s3) * (s0 - s4) * (s1 - s4) +
                         s2 * s2 * (s2 - s3) *
                             (-B * s0 * s0 * (s0 - s2) * (s0 - s3) * (s0 - s4) +
                              A * s1 * s1 * (s1 - s2) * (s1 - s3) * (s1 - s4)) *
                             (s2 - s4)) *
                  (s3 - s4)) *
                 s4 * s4) /
            ((s0 - s1) * (s0 - s2) * (s1 - s2) * (s0 - s3) * (s1 - s3) *
             (s2 - s3) * (s0 - s4) * (s1 - s4) * (s2 - s4) * (s3 - s4) *
             (s1 * s2 * s3 * s4 +
              s0 * (s2 * s3 * s4 + s1 * (s3 * s4 + s2 * (s3 + s4)))));

  c[1].re =
      (-C * s0 * s0 * s0 * s0 * s0 * s1 * s1 * s1 * s1 * s3 * s3 * s3 +
       C * s0 * s0 * s0 * s0 * s1 * s1 * s1 * s1 * s1 * s3 * s3 * s3 +
       B * s0 * s0 * s0 * s0 * s0 * s2 * s2 * s2 * s2 * s3 * s3 * s3 -
       A * s1 * s1 * s1 * s1 * s1 * s2 * s2 * s2 * s2 * s3 * s3 * s3 -
       B * s0 * s0 * s0 * s0 * s2 * s2 * s2 * s2 * s2 * s3 * s3 * s3 +
       A * s1 * s1 * s1 * s1 * s2 * s2 * s2 * s2 * s2 * s3 * s3 * s3 +
       C * s0 * s0 * s0 * s0 * s0 * s1 * s1 * s1 * s3 * s3 * s3 * s3 -
       C * s0 * s0 * s0 * s1 * s1 * s1 * s1 * s1 * s3 * s3 * s3 * s3 -
       B * s0 * s0 * s0 * s0 * s0 * s2 * s2 * s2 * s3 * s3 * s3 * s3 +
       A * s1 * s1 * s1 * s1 * s1 * s2 * s2 * s2 * s3 * s3 * s3 * s3 +
       B * s0 * s0 * s0 * s2 * s2 * s2 * s2 * s2 * s3 * s3 * s3 * s3 -
       A * s1 * s1 * s1 * s2 * s2 * s2 * s2 * s2 * s3 * s3 * s3 * s3 -
       C * s0 * s0 * s0 * s0 * s1 * s1 * s1 * s3 * s3 * s3 * s3 * s3 +
       C * s0 * s0 * s0 * s1 * s1 * s1 * s1 * s3 * s3 * s3 * s3 * s3 +
       B * s0 * s0 * s0 * s0 * s2 * s2 * s2 * s3 * s3 * s3 * s3 * s3 -
       A * s1 * s1 * s1 * s1 * s2 * s2 * s2 * s3 * s3 * s3 * s3 * s3 -
       B * s0 * s0 * s0 * s2 * s2 * s2 * s2 * s3 * s3 * s3 * s3 * s3 +
       A * s1 * s1 * s1 * s2 * s2 * s2 * s2 * s3 * s3 * s3 * s3 * s3 +
       E * (s1 * s1 * s1 * (s1 - s2) * s2 * s2 * s2 * (s1 - s3) * (s2 - s3) *
                s3 * s3 * s3 +
            s0 * s0 * s0 * s0 * s0 *
                (s2 * s2 * s2 * s3 * s3 * s3 * (-s2 + s3) +
                 s1 * s1 * s1 * s1 * (-s2 * s2 * s2 + s3 * s3 * s3) +
                 s1 * s1 * s1 * (s2 * s2 * s2 * s2 - s3 * s3 * s3 * s3)) +
            s0 * s0 * s0 * (s2 * s2 * s2 * s2 * s3 * s3 * s3 * s3 * (-s2 + s3) +
                            s1 * s1 * s1 * s1 * s1 *
                                (-s2 * s2 * s2 * s2 + s3 * s3 * s3 * s3) +
                            s1 * s1 * s1 * s1 * (s2 * s2 * s2 * s2 * s2 -
                                                 s3 * s3 * s3 * s3 * s3)) +
            s0 * s0 * s0 * s0 *
                (s2 * s2 * s2 * s3 * s3 * s3 * (s2 * s2 - s3 * s3) +
                 s1 * s1 * s1 * s1 * s1 * (s2 * s2 * s2 - s3 * s3 * s3) +
                 s1 * s1 * s1 *
                     (-s2 * s2 * s2 * s2 * s2 + s3 * s3 * s3 * s3 * s3))) +
       C * s0 * s0 * s0 * s0 * s0 * s1 * s1 * s1 * s1 * s4 * s4 * s4 -
       C * s0 * s0 * s0 * s0 * s1 * s1 * s1 * s1 * s1 * s4 * s4 * s4 -
       B * s0 * s0 * s0 * s0 * s0 * s2 * s2 * s2 * s2 * s4 * s4 * s4 +
       A * s1 * s1 * s1 * s1 * s1 * s2 * s2 * s2 * s2 * s4 * s4 * s4 +
       B * s0 * s0 * s0 * s0 * s2 * s2 * s2 * s2 * s2 * s4 * s4 * s4 -
       A * s1 * s1 * s1 * s1 * s2 * s2 * s2 * s2 * s2 * s4 * s4 * s4 +
       B * s0 * s0 * s0 * s0 * s0 * s3 * s3 * s3 * s3 * s4 * s4 * s4 -
       C * s0 * s0 * s0 * s0 * s0 * s3 * s3 * s3 * s3 * s4 * s4 * s4 -
       A * s1 * s1 * s1 * s1 * s1 * s3 * s3 * s3 * s3 * s4 * s4 * s4 +
       C * s1 * s1 * s1 * s1 * s1 * s3 * s3 * s3 * s3 * s4 * s4 * s4 +
       A * s2 * s2 * s2 * s2 * s2 * s3 * s3 * s3 * s3 * s4 * s4 * s4 -
       B * s2 * s2 * s2 * s2 * s2 * s3 * s3 * s3 * s3 * s4 * s4 * s4 -
       B * s0 * s0 * s0 * s0 * s3 * s3 * s3 * s3 * s3 * s4 * s4 * s4 +
       C * s0 * s0 * s0 * s0 * s3 * s3 * s3 * s3 * s3 * s4 * s4 * s4 +
       A * s1 * s1 * s1 * s1 * s3 * s3 * s3 * s3 * s3 * s4 * s4 * s4 -
       C * s1 * s1 * s1 * s1 * s3 * s3 * s3 * s3 * s3 * s4 * s4 * s4 -
       A * s2 * s2 * s2 * s2 * s3 * s3 * s3 * s3 * s3 * s4 * s4 * s4 +
       B * s2 * s2 * s2 * s2 * s3 * s3 * s3 * s3 * s3 * s4 * s4 * s4 -
       C * s0 * s0 * s0 * s0 * s0 * s1 * s1 * s1 * s4 * s4 * s4 * s4 +
       C * s0 * s0 * s0 * s1 * s1 * s1 * s1 * s1 * s4 * s4 * s4 * s4 +
       B * s0 * s0 * s0 * s0 * s0 * s2 * s2 * s2 * s4 * s4 * s4 * s4 -
       A * s1 * s1 * s1 * s1 * s1 * s2 * s2 * s2 * s4 * s4 * s4 * s4 -
       B * s0 * s0 * s0 * s2 * s2 * s2 * s2 * s2 * s4 * s4 * s4 * s4 +
       A * s1 * s1 * s1 * s2 * s2 * s2 * s2 * s2 * s4 * s4 * s4 * s4 -
       B * s0 * s0 * s0 * s0 * s0 * s3 * s3 * s3 * s4 * s4 * s4 * s4 +
       C * s0 * s0 * s0 * s0 * s0 * s3 * s3 * s3 * s4 * s4 * s4 * s4 +
       A * s1 * s1 * s1 * s1 * s1 * s3 * s3 * s3 * s4 * s4 * s4 * s4 -
       C * s1 * s1 * s1 * s1 * s1 * s3 * s3 * s3 * s4 * s4 * s4 * s4 -
       A * s2 * s2 * s2 * s2 * s2 * s3 * s3 * s3 * s4 * s4 * s4 * s4 +
       B * s2 * s2 * s2 * s2 * s2 * s3 * s3 * s3 * s4 * s4 * s4 * s4 +
       B * s0 * s0 * s0 * s3 * s3 * s3 * s3 * s3 * s4 * s4 * s4 * s4 -
       C * s0 * s0 * s0 * s3 * s3 * s3 * s3 * s3 * s4 * s4 * s4 * s4 -
       A * s1 * s1 * s1 * s3 * s3 * s3 * s3 * s3 * s4 * s4 * s4 * s4 +
       C * s1 * s1 * s1 * s3 * s3 * s3 * s3 * s3 * s4 * s4 * s4 * s4 +
       A * s2 * s2 * s2 * s3 * s3 * s3 * s3 * s3 * s4 * s4 * s4 * s4 -
       B * s2 * s2 * s2 * s3 * s3 * s3 * s3 * s3 * s4 * s4 * s4 * s4 +
       C * s0 * s0 * s0 * s0 * s1 * s1 * s1 * s4 * s4 * s4 * s4 * s4 -
       C * s0 * s0 * s0 * s1 * s1 * s1 * s1 * s4 * s4 * s4 * s4 * s4 -
       B * s0 * s0 * s0 * s0 * s2 * s2 * s2 * s4 * s4 * s4 * s4 * s4 +
       A * s1 * s1 * s1 * s1 * s2 * s2 * s2 * s4 * s4 * s4 * s4 * s4 +
       B * s0 * s0 * s0 * s2 * s2 * s2 * s2 * s4 * s4 * s4 * s4 * s4 -
       A * s1 * s1 * s1 * s2 * s2 * s2 * s2 * s4 * s4 * s4 * s4 * s4 +
       B * s0 * s0 * s0 * s0 * s3 * s3 * s3 * s4 * s4 * s4 * s4 * s4 -
       C * s0 * s0 * s0 * s0 * s3 * s3 * s3 * s4 * s4 * s4 * s4 * s4 -
       A * s1 * s1 * s1 * s1 * s3 * s3 * s3 * s4 * s4 * s4 * s4 * s4 +
       C * s1 * s1 * s1 * s1 * s3 * s3 * s3 * s4 * s4 * s4 * s4 * s4 +
       A * s2 * s2 * s2 * s2 * s3 * s3 * s3 * s4 * s4 * s4 * s4 * s4 -
       B * s2 * s2 * s2 * s2 * s3 * s3 * s3 * s4 * s4 * s4 * s4 * s4 -
       B * s0 * s0 * s0 * s3 * s3 * s3 * s3 * s4 * s4 * s4 * s4 * s4 +
       C * s0 * s0 * s0 * s3 * s3 * s3 * s3 * s4 * s4 * s4 * s4 * s4 +
       A * s1 * s1 * s1 * s3 * s3 * s3 * s3 * s4 * s4 * s4 * s4 * s4 -
       C * s1 * s1 * s1 * s3 * s3 * s3 * s3 * s4 * s4 * s4 * s4 * s4 -
       A * s2 * s2 * s2 * s3 * s3 * s3 * s3 * s4 * s4 * s4 * s4 * s4 +
       B * s2 * s2 * s2 * s3 * s3 * s3 * s3 * s4 * s4 * s4 * s4 * s4 +
       D * (-s1 * s1 * s1 * (s1 - s2) * s2 * s2 * s2 * (s1 - s4) * (s2 - s4) *
                s4 * s4 * s4 +
            s0 * s0 * s0 * s0 * s0 *
                (s2 * s2 * s2 * (s2 - s4) * s4 * s4 * s4 +
                 s1 * s1 * s1 * s1 * (s2 * s2 * s2 - s4 * s4 * s4) +
                 s1 * s1 * s1 * (-s2 * s2 * s2 * s2 + s4 * s4 * s4 * s4)) +
            s0 * s0 * s0 * s0 *
                (-s2 * s2 * s2 * s2 * s2 * s4 * s4 * s4 +
                 s2 * s2 * s2 * s4 * s4 * s4 * s4 * s4 +
                 s1 * s1 * s1 * s1 * s1 * (-s2 * s2 * s2 + s4 * s4 * s4) +
                 s1 * s1 * s1 *
                     (s2 * s2 * s2 * s2 * s2 - s4 * s4 * s4 * s4 * s4)) +
            s0 * s0 * s0 * (s2 * s2 * s2 * s2 * (s2 - s4) * s4 * s4 * s4 * s4 +
                            s1 * s1 * s1 * s1 * s1 *
                                (s2 * s2 * s2 * s2 - s4 * s4 * s4 * s4) +
                            s1 * s1 * s1 * s1 * (-s2 * s2 * s2 * s2 * s2 +
                                                 s4 * s4 * s4 * s4 * s4)))) /
      ((s0 - s1) * (s0 - s2) * (s1 - s2) * (s0 - s3) * (s1 - s3) * (s2 - s3) *
       (s0 - s4) * (s1 - s4) * (s2 - s4) * (s3 - s4) *
       (s1 * s2 * s3 * s4 +
        s0 * (s2 * s3 * s4 + s1 * (s3 * s4 + s2 * (s3 + s4)))));

  c[2].re =
      (C * s0 * s0 * s0 * s0 * s0 * s1 * s1 * s1 * s1 * s3 * s3 -
       C * s0 * s0 * s0 * s0 * s1 * s1 * s1 * s1 * s1 * s3 * s3 -
       B * s0 * s0 * s0 * s0 * s0 * s2 * s2 * s2 * s2 * s3 * s3 +
       A * s1 * s1 * s1 * s1 * s1 * s2 * s2 * s2 * s2 * s3 * s3 +
       B * s0 * s0 * s0 * s0 * s2 * s2 * s2 * s2 * s2 * s3 * s3 -
       A * s1 * s1 * s1 * s1 * s2 * s2 * s2 * s2 * s2 * s3 * s3 -
       C * s0 * s0 * s0 * s0 * s0 * s1 * s1 * s3 * s3 * s3 * s3 +
       C * s0 * s0 * s1 * s1 * s1 * s1 * s1 * s3 * s3 * s3 * s3 +
       B * s0 * s0 * s0 * s0 * s0 * s2 * s2 * s3 * s3 * s3 * s3 -
       A * s1 * s1 * s1 * s1 * s1 * s2 * s2 * s3 * s3 * s3 * s3 -
       B * s0 * s0 * s2 * s2 * s2 * s2 * s2 * s3 * s3 * s3 * s3 +
       A * s1 * s1 * s2 * s2 * s2 * s2 * s2 * s3 * s3 * s3 * s3 +
       C * s0 * s0 * s0 * s0 * s1 * s1 * s3 * s3 * s3 * s3 * s3 -
       C * s0 * s0 * s1 * s1 * s1 * s1 * s3 * s3 * s3 * s3 * s3 -
       B * s0 * s0 * s0 * s0 * s2 * s2 * s3 * s3 * s3 * s3 * s3 +
       A * s1 * s1 * s1 * s1 * s2 * s2 * s3 * s3 * s3 * s3 * s3 +
       B * s0 * s0 * s2 * s2 * s2 * s2 * s3 * s3 * s3 * s3 * s3 -
       A * s1 * s1 * s2 * s2 * s2 * s2 * s3 * s3 * s3 * s3 * s3 +
       E * (s0 * s0 * s0 * s0 * s0 * (s1 * s1 - s2 * s2) * (s1 * s1 - s3 * s3) *
                (s2 * s2 - s3 * s3) +
            s1 * s1 * s2 * s2 * s3 * s3 *
                (s2 * s2 * s3 * s3 * (-s2 + s3) +
                 s1 * s1 * s1 * (-s2 * s2 + s3 * s3) +
                 s1 * s1 * (s2 * s2 * s2 - s3 * s3 * s3)) +
            s0 * s0 * s0 * s0 *
                (-s2 * s2 * s2 * s2 * s2 * s3 * s3 +
                 s2 * s2 * s3 * s3 * s3 * s3 * s3 +
                 s1 * s1 * s1 * s1 * s1 * (-s2 * s2 + s3 * s3) +
                 s1 * s1 * (s2 * s2 * s2 * s2 * s2 - s3 * s3 * s3 * s3 * s3)) +
            s0 * s0 * (s2 * s2 * s2 * s2 * (s2 - s3) * s3 * s3 * s3 * s3 +
                       s1 * s1 * s1 * s1 * s1 *
                           (s2 * s2 * s2 * s2 - s3 * s3 * s3 * s3) +
                       s1 * s1 * s1 * s1 * (-s2 * s2 * s2 * s2 * s2 +
                                            s3 * s3 * s3 * s3 * s3))) -
       C * s0 * s0 * s0 * s0 * s0 * s1 * s1 * s1 * s1 * s4 * s4 +
       C * s0 * s0 * s0 * s0 * s1 * s1 * s1 * s1 * s1 * s4 * s4 +
       B * s0 * s0 * s0 * s0 * s0 * s2 * s2 * s2 * s2 * s4 * s4 -
       A * s1 * s1 * s1 * s1 * s1 * s2 * s2 * s2 * s2 * s4 * s4 -
       B * s0 * s0 * s0 * s0 * s2 * s2 * s2 * s2 * s2 * s4 * s4 +
       A * s1 * s1 * s1 * s1 * s2 * s2 * s2 * s2 * s2 * s4 * s4 -
       B * s0 * s0 * s0 * s0 * s0 * s3 * s3 * s3 * s3 * s4 * s4 +
       C * s0 * s0 * s0 * s0 * s0 * s3 * s3 * s3 * s3 * s4 * s4 +
       A * s1 * s1 * s1 * s1 * s1 * s3 * s3 * s3 * s3 * s4 * s4 -
       C * s1 * s1 * s1 * s1 * s1 * s3 * s3 * s3 * s3 * s4 * s4 -
       A * s2 * s2 * s2 * s2 * s2 * s3 * s3 * s3 * s3 * s4 * s4 +
       B * s2 * s2 * s2 * s2 * s2 * s3 * s3 * s3 * s3 * s4 * s4 +
       B * s0 * s0 * s0 * s0 * s3 * s3 * s3 * s3 * s3 * s4 * s4 -
       C * s0 * s0 * s0 * s0 * s3 * s3 * s3 * s3 * s3 * s4 * s4 -
       A * s1 * s1 * s1 * s1 * s3 * s3 * s3 * s3 * s3 * s4 * s4 +
       C * s1 * s1 * s1 * s1 * s3 * s3 * s3 * s3 * s3 * s4 * s4 +
       A * s2 * s2 * s2 * s2 * s3 * s3 * s3 * s3 * s3 * s4 * s4 -
       B * s2 * s2 * s2 * s2 * s3 * s3 * s3 * s3 * s3 * s4 * s4 +
       C * s0 * s0 * s0 * s0 * s0 * s1 * s1 * s4 * s4 * s4 * s4 -
       C * s0 * s0 * s1 * s1 * s1 * s1 * s1 * s4 * s4 * s4 * s4 -
       B * s0 * s0 * s0 * s0 * s0 * s2 * s2 * s4 * s4 * s4 * s4 +
       A * s1 * s1 * s1 * s1 * s1 * s2 * s2 * s4 * s4 * s4 * s4 +
       B * s0 * s0 * s2 * s2 * s2 * s2 * s2 * s4 * s4 * s4 * s4 -
       A * s1 * s1 * s2 * s2 * s2 * s2 * s2 * s4 * s4 * s4 * s4 +
       B * s0 * s0 * s0 * s0 * s0 * s3 * s3 * s4 * s4 * s4 * s4 -
       C * s0 * s0 * s0 * s0 * s0 * s3 * s3 * s4 * s4 * s4 * s4 -
       A * s1 * s1 * s1 * s1 * s1 * s3 * s3 * s4 * s4 * s4 * s4 +
       C * s1 * s1 * s1 * s1 * s1 * s3 * s3 * s4 * s4 * s4 * s4 +
       A * s2 * s2 * s2 * s2 * s2 * s3 * s3 * s4 * s4 * s4 * s4 -
       B * s2 * s2 * s2 * s2 * s2 * s3 * s3 * s4 * s4 * s4 * s4 -
       B * s0 * s0 * s3 * s3 * s3 * s3 * s3 * s4 * s4 * s4 * s4 +
       C * s0 * s0 * s3 * s3 * s3 * s3 * s3 * s4 * s4 * s4 * s4 +
       A * s1 * s1 * s3 * s3 * s3 * s3 * s3 * s4 * s4 * s4 * s4 -
       C * s1 * s1 * s3 * s3 * s3 * s3 * s3 * s4 * s4 * s4 * s4 -
       A * s2 * s2 * s3 * s3 * s3 * s3 * s3 * s4 * s4 * s4 * s4 +
       B * s2 * s2 * s3 * s3 * s3 * s3 * s3 * s4 * s4 * s4 * s4 -
       C * s0 * s0 * s0 * s0 * s1 * s1 * s4 * s4 * s4 * s4 * s4 +
       C * s0 * s0 * s1 * s1 * s1 * s1 * s4 * s4 * s4 * s4 * s4 +
       B * s0 * s0 * s0 * s0 * s2 * s2 * s4 * s4 * s4 * s4 * s4 -
       A * s1 * s1 * s1 * s1 * s2 * s2 * s4 * s4 * s4 * s4 * s4 -
       B * s0 * s0 * s2 * s2 * s2 * s2 * s4 * s4 * s4 * s4 * s4 +
       A * s1 * s1 * s2 * s2 * s2 * s2 * s4 * s4 * s4 * s4 * s4 -
       B * s0 * s0 * s0 * s0 * s3 * s3 * s4 * s4 * s4 * s4 * s4 +
       C * s0 * s0 * s0 * s0 * s3 * s3 * s4 * s4 * s4 * s4 * s4 +
       A * s1 * s1 * s1 * s1 * s3 * s3 * s4 * s4 * s4 * s4 * s4 -
       C * s1 * s1 * s1 * s1 * s3 * s3 * s4 * s4 * s4 * s4 * s4 -
       A * s2 * s2 * s2 * s2 * s3 * s3 * s4 * s4 * s4 * s4 * s4 +
       B * s2 * s2 * s2 * s2 * s3 * s3 * s4 * s4 * s4 * s4 * s4 +
       B * s0 * s0 * s3 * s3 * s3 * s3 * s4 * s4 * s4 * s4 * s4 -
       C * s0 * s0 * s3 * s3 * s3 * s3 * s4 * s4 * s4 * s4 * s4 -
       A * s1 * s1 * s3 * s3 * s3 * s3 * s4 * s4 * s4 * s4 * s4 +
       C * s1 * s1 * s3 * s3 * s3 * s3 * s4 * s4 * s4 * s4 * s4 +
       A * s2 * s2 * s3 * s3 * s3 * s3 * s4 * s4 * s4 * s4 * s4 -
       B * s2 * s2 * s3 * s3 * s3 * s3 * s4 * s4 * s4 * s4 * s4 +
       D * (-s0 * s0 * s0 * s0 * s0 * (s1 * s1 - s2 * s2) *
                (s1 * s1 - s4 * s4) * (s2 * s2 - s4 * s4) +
            s1 * s1 * s2 * s2 * s4 * s4 *
                (s2 * s2 * (s2 - s4) * s4 * s4 +
                 s1 * s1 * s1 * (s2 * s2 - s4 * s4) +
                 s1 * s1 * (-s2 * s2 * s2 + s4 * s4 * s4)) +
            s0 * s0 * (s2 * s2 * s2 * s2 * s4 * s4 * s4 * s4 * (-s2 + s4) +
                       s1 * s1 * s1 * s1 * s1 *
                           (-s2 * s2 * s2 * s2 + s4 * s4 * s4 * s4) +
                       s1 * s1 * s1 * s1 *
                           (s2 * s2 * s2 * s2 * s2 - s4 * s4 * s4 * s4 * s4)) +
            s0 * s0 * s0 * s0 *
                (s1 * s1 * s1 * s1 * s1 * (s2 * s2 - s4 * s4) +
                 s2 * s2 * s4 * s4 * (s2 * s2 * s2 - s4 * s4 * s4) +
                 s1 * s1 *
                     (-s2 * s2 * s2 * s2 * s2 + s4 * s4 * s4 * s4 * s4)))) /
      ((s0 - s1) * (s0 - s2) * (s1 - s2) * (s0 - s3) * (s1 - s3) * (s2 - s3) *
       (s0 - s4) * (s1 - s4) * (s2 - s4) * (s3 - s4) *
       (s1 * s2 * s3 * s4 +
        s0 * (s2 * s3 * s4 + s1 * (s3 * s4 + s2 * (s3 + s4)))));

  c[3].re =
      (-C * s0 * s0 * s0 * s0 * s0 * s1 * s1 * s1 * s3 * s3 +
       C * s0 * s0 * s0 * s1 * s1 * s1 * s1 * s1 * s3 * s3 +
       B * s0 * s0 * s0 * s0 * s0 * s2 * s2 * s2 * s3 * s3 -
       A * s1 * s1 * s1 * s1 * s1 * s2 * s2 * s2 * s3 * s3 -
       B * s0 * s0 * s0 * s2 * s2 * s2 * s2 * s2 * s3 * s3 +
       A * s1 * s1 * s1 * s2 * s2 * s2 * s2 * s2 * s3 * s3 +
       C * s0 * s0 * s0 * s0 * s0 * s1 * s1 * s3 * s3 * s3 -
       C * s0 * s0 * s1 * s1 * s1 * s1 * s1 * s3 * s3 * s3 -
       B * s0 * s0 * s0 * s0 * s0 * s2 * s2 * s3 * s3 * s3 +
       A * s1 * s1 * s1 * s1 * s1 * s2 * s2 * s3 * s3 * s3 +
       B * s0 * s0 * s2 * s2 * s2 * s2 * s2 * s3 * s3 * s3 -
       A * s1 * s1 * s2 * s2 * s2 * s2 * s2 * s3 * s3 * s3 -
       C * s0 * s0 * s0 * s1 * s1 * s3 * s3 * s3 * s3 * s3 +
       C * s0 * s0 * s1 * s1 * s1 * s3 * s3 * s3 * s3 * s3 +
       B * s0 * s0 * s0 * s2 * s2 * s3 * s3 * s3 * s3 * s3 -
       A * s1 * s1 * s1 * s2 * s2 * s3 * s3 * s3 * s3 * s3 -
       B * s0 * s0 * s2 * s2 * s2 * s3 * s3 * s3 * s3 * s3 +
       A * s1 * s1 * s2 * s2 * s2 * s3 * s3 * s3 * s3 * s3 +
       E * (s0 * s0 * s0 * s0 * s0 * (s2 * s2 * s3 * s3 * (-s2 + s3) +
                                      s1 * s1 * s1 * (-s2 * s2 + s3 * s3) +
                                      s1 * s1 * (s2 * s2 * s2 - s3 * s3 * s3)) +
            s1 * s1 * s2 * s2 * s3 * s3 *
                (s1 * s1 * s1 * (s2 - s3) + s2 * s3 * (s2 * s2 - s3 * s3) +
                 s1 * (-s2 * s2 * s2 + s3 * s3 * s3)) +
            s0 * s0 * (-s2 * s2 * s2 * s2 * s2 * s3 * s3 * s3 +
                       s2 * s2 * s2 * s3 * s3 * s3 * s3 * s3 +
                       s1 * s1 * s1 * s1 * s1 * (-s2 * s2 * s2 + s3 * s3 * s3) +
                       s1 * s1 * s1 *
                           (s2 * s2 * s2 * s2 * s2 - s3 * s3 * s3 * s3 * s3)) +
            s0 * s0 * s0 * (s1 * s1 * s1 * s1 * s1 * (s2 * s2 - s3 * s3) +
                            s2 * s2 * s3 * s3 * (s2 * s2 * s2 - s3 * s3 * s3) +
                            s1 * s1 * (-s2 * s2 * s2 * s2 * s2 +
                                       s3 * s3 * s3 * s3 * s3))) +
       C * s0 * s0 * s0 * s0 * s0 * s1 * s1 * s1 * s4 * s4 -
       C * s0 * s0 * s0 * s1 * s1 * s1 * s1 * s1 * s4 * s4 -
       B * s0 * s0 * s0 * s0 * s0 * s2 * s2 * s2 * s4 * s4 +
       A * s1 * s1 * s1 * s1 * s1 * s2 * s2 * s2 * s4 * s4 +
       B * s0 * s0 * s0 * s2 * s2 * s2 * s2 * s2 * s4 * s4 -
       A * s1 * s1 * s1 * s2 * s2 * s2 * s2 * s2 * s4 * s4 +
       B * s0 * s0 * s0 * s0 * s0 * s3 * s3 * s3 * s4 * s4 -
       C * s0 * s0 * s0 * s0 * s0 * s3 * s3 * s3 * s4 * s4 -
       A * s1 * s1 * s1 * s1 * s1 * s3 * s3 * s3 * s4 * s4 +
       C * s1 * s1 * s1 * s1 * s1 * s3 * s3 * s3 * s4 * s4 +
       A * s2 * s2 * s2 * s2 * s2 * s3 * s3 * s3 * s4 * s4 -
       B * s2 * s2 * s2 * s2 * s2 * s3 * s3 * s3 * s4 * s4 -
       B * s0 * s0 * s0 * s3 * s3 * s3 * s3 * s3 * s4 * s4 +
       C * s0 * s0 * s0 * s3 * s3 * s3 * s3 * s3 * s4 * s4 +
       A * s1 * s1 * s1 * s3 * s3 * s3 * s3 * s3 * s4 * s4 -
       C * s1 * s1 * s1 * s3 * s3 * s3 * s3 * s3 * s4 * s4 -
       A * s2 * s2 * s2 * s3 * s3 * s3 * s3 * s3 * s4 * s4 +
       B * s2 * s2 * s2 * s3 * s3 * s3 * s3 * s3 * s4 * s4 -
       C * s0 * s0 * s0 * s0 * s0 * s1 * s1 * s4 * s4 * s4 +
       C * s0 * s0 * s1 * s1 * s1 * s1 * s1 * s4 * s4 * s4 +
       B * s0 * s0 * s0 * s0 * s0 * s2 * s2 * s4 * s4 * s4 -
       A * s1 * s1 * s1 * s1 * s1 * s2 * s2 * s4 * s4 * s4 -
       B * s0 * s0 * s2 * s2 * s2 * s2 * s2 * s4 * s4 * s4 +
       A * s1 * s1 * s2 * s2 * s2 * s2 * s2 * s4 * s4 * s4 -
       B * s0 * s0 * s0 * s0 * s0 * s3 * s3 * s4 * s4 * s4 +
       C * s0 * s0 * s0 * s0 * s0 * s3 * s3 * s4 * s4 * s4 +
       A * s1 * s1 * s1 * s1 * s1 * s3 * s3 * s4 * s4 * s4 -
       C * s1 * s1 * s1 * s1 * s1 * s3 * s3 * s4 * s4 * s4 -
       A * s2 * s2 * s2 * s2 * s2 * s3 * s3 * s4 * s4 * s4 +
       B * s2 * s2 * s2 * s2 * s2 * s3 * s3 * s4 * s4 * s4 +
       B * s0 * s0 * s3 * s3 * s3 * s3 * s3 * s4 * s4 * s4 -
       C * s0 * s0 * s3 * s3 * s3 * s3 * s3 * s4 * s4 * s4 -
       A * s1 * s1 * s3 * s3 * s3 * s3 * s3 * s4 * s4 * s4 +
       C * s1 * s1 * s3 * s3 * s3 * s3 * s3 * s4 * s4 * s4 +
       A * s2 * s2 * s3 * s3 * s3 * s3 * s3 * s4 * s4 * s4 -
       B * s2 * s2 * s3 * s3 * s3 * s3 * s3 * s4 * s4 * s4 +
       C * s0 * s0 * s0 * s1 * s1 * s4 * s4 * s4 * s4 * s4 -
       C * s0 * s0 * s1 * s1 * s1 * s4 * s4 * s4 * s4 * s4 -
       B * s0 * s0 * s0 * s2 * s2 * s4 * s4 * s4 * s4 * s4 +
       A * s1 * s1 * s1 * s2 * s2 * s4 * s4 * s4 * s4 * s4 +
       B * s0 * s0 * s2 * s2 * s2 * s4 * s4 * s4 * s4 * s4 -
       A * s1 * s1 * s2 * s2 * s2 * s4 * s4 * s4 * s4 * s4 +
       B * s0 * s0 * s0 * s3 * s3 * s4 * s4 * s4 * s4 * s4 -
       C * s0 * s0 * s0 * s3 * s3 * s4 * s4 * s4 * s4 * s4 -
       A * s1 * s1 * s1 * s3 * s3 * s4 * s4 * s4 * s4 * s4 +
       C * s1 * s1 * s1 * s3 * s3 * s4 * s4 * s4 * s4 * s4 +
       A * s2 * s2 * s2 * s3 * s3 * s4 * s4 * s4 * s4 * s4 -
       B * s2 * s2 * s2 * s3 * s3 * s4 * s4 * s4 * s4 * s4 -
       B * s0 * s0 * s3 * s3 * s3 * s4 * s4 * s4 * s4 * s4 +
       C * s0 * s0 * s3 * s3 * s3 * s4 * s4 * s4 * s4 * s4 +
       A * s1 * s1 * s3 * s3 * s3 * s4 * s4 * s4 * s4 * s4 -
       C * s1 * s1 * s3 * s3 * s3 * s4 * s4 * s4 * s4 * s4 -
       A * s2 * s2 * s3 * s3 * s3 * s4 * s4 * s4 * s4 * s4 +
       B * s2 * s2 * s3 * s3 * s3 * s4 * s4 * s4 * s4 * s4 +
       D * (s1 * s1 * s2 * s2 * s4 * s4 *
                (-s2 * s2 * s2 * s4 + s2 * s4 * s4 * s4 +
                 s1 * s1 * s1 * (-s2 + s4) +
                 s1 * (s2 * s2 * s2 - s4 * s4 * s4)) +
            s0 * s0 * s0 * s0 * s0 *
                (s2 * s2 * (s2 - s4) * s4 * s4 +
                 s1 * s1 * s1 * (s2 * s2 - s4 * s4) +
                 s1 * s1 * (-s2 * s2 * s2 + s4 * s4 * s4)) +
            s0 * s0 * s0 *
                (-s2 * s2 * s2 * s2 * s2 * s4 * s4 +
                 s2 * s2 * s4 * s4 * s4 * s4 * s4 +
                 s1 * s1 * s1 * s1 * s1 * (-s2 * s2 + s4 * s4) +
                 s1 * s1 * (s2 * s2 * s2 * s2 * s2 - s4 * s4 * s4 * s4 * s4)) +
            s0 * s0 * (s2 * s2 * s2 * s4 * s4 * s4 * (s2 * s2 - s4 * s4) +
                       s1 * s1 * s1 * s1 * s1 * (s2 * s2 * s2 - s4 * s4 * s4) +
                       s1 * s1 * s1 * (-s2 * s2 * s2 * s2 * s2 +
                                       s4 * s4 * s4 * s4 * s4)))) /
      ((s0 - s1) * (s0 - s2) * (s1 - s2) * (s0 - s3) * (s1 - s3) * (s2 - s3) *
       (s0 - s4) * (s1 - s4) * (s2 - s4) * (s3 - s4) *
       (s1 * s2 * s3 * s4 +
        s0 * (s2 * s3 * s4 + s1 * (s3 * s4 + s2 * (s3 + s4)))));

  c[4].re = (C * s0 * s0 * s0 * s0 * s1 * s1 * s1 * s3 * s3 -
             C * s0 * s0 * s0 * s1 * s1 * s1 * s1 * s3 * s3 -
             B * s0 * s0 * s0 * s0 * s2 * s2 * s2 * s3 * s3 +
             A * s1 * s1 * s1 * s1 * s2 * s2 * s2 * s3 * s3 +
             B * s0 * s0 * s0 * s2 * s2 * s2 * s2 * s3 * s3 -
             A * s1 * s1 * s1 * s2 * s2 * s2 * s2 * s3 * s3 -
             C * s0 * s0 * s0 * s0 * s1 * s1 * s3 * s3 * s3 +
             C * s0 * s0 * s1 * s1 * s1 * s1 * s3 * s3 * s3 +
             B * s0 * s0 * s0 * s0 * s2 * s2 * s3 * s3 * s3 -
             A * s1 * s1 * s1 * s1 * s2 * s2 * s3 * s3 * s3 -
             B * s0 * s0 * s2 * s2 * s2 * s2 * s3 * s3 * s3 +
             A * s1 * s1 * s2 * s2 * s2 * s2 * s3 * s3 * s3 +
             C * s0 * s0 * s0 * s1 * s1 * s3 * s3 * s3 * s3 -
             C * s0 * s0 * s1 * s1 * s1 * s3 * s3 * s3 * s3 -
             B * s0 * s0 * s0 * s2 * s2 * s3 * s3 * s3 * s3 +
             A * s1 * s1 * s1 * s2 * s2 * s3 * s3 * s3 * s3 +
             B * s0 * s0 * s2 * s2 * s2 * s3 * s3 * s3 * s3 -
             A * s1 * s1 * s2 * s2 * s2 * s3 * s3 * s3 * s3 +
             E * (s0 - s1) * (s0 - s2) * (s1 - s2) * (s0 - s3) * (s1 - s3) *
                 (s2 - s3) * (s1 * s2 * s3 + s0 * (s2 * s3 + s1 * (s2 + s3))) -
             C * s0 * s0 * s0 * s0 * s1 * s1 * s1 * s4 * s4 +
             C * s0 * s0 * s0 * s1 * s1 * s1 * s1 * s4 * s4 +
             B * s0 * s0 * s0 * s0 * s2 * s2 * s2 * s4 * s4 -
             A * s1 * s1 * s1 * s1 * s2 * s2 * s2 * s4 * s4 -
             B * s0 * s0 * s0 * s2 * s2 * s2 * s2 * s4 * s4 +
             A * s1 * s1 * s1 * s2 * s2 * s2 * s2 * s4 * s4 -
             B * s0 * s0 * s0 * s0 * s3 * s3 * s3 * s4 * s4 +
             C * s0 * s0 * s0 * s0 * s3 * s3 * s3 * s4 * s4 +
             A * s1 * s1 * s1 * s1 * s3 * s3 * s3 * s4 * s4 -
             C * s1 * s1 * s1 * s1 * s3 * s3 * s3 * s4 * s4 -
             A * s2 * s2 * s2 * s2 * s3 * s3 * s3 * s4 * s4 +
             B * s2 * s2 * s2 * s2 * s3 * s3 * s3 * s4 * s4 +
             B * s0 * s0 * s0 * s3 * s3 * s3 * s3 * s4 * s4 -
             C * s0 * s0 * s0 * s3 * s3 * s3 * s3 * s4 * s4 -
             A * s1 * s1 * s1 * s3 * s3 * s3 * s3 * s4 * s4 +
             C * s1 * s1 * s1 * s3 * s3 * s3 * s3 * s4 * s4 +
             A * s2 * s2 * s2 * s3 * s3 * s3 * s3 * s4 * s4 -
             B * s2 * s2 * s2 * s3 * s3 * s3 * s3 * s4 * s4 +
             C * s0 * s0 * s0 * s0 * s1 * s1 * s4 * s4 * s4 -
             C * s0 * s0 * s1 * s1 * s1 * s1 * s4 * s4 * s4 -
             B * s0 * s0 * s0 * s0 * s2 * s2 * s4 * s4 * s4 +
             A * s1 * s1 * s1 * s1 * s2 * s2 * s4 * s4 * s4 +
             B * s0 * s0 * s2 * s2 * s2 * s2 * s4 * s4 * s4 -
             A * s1 * s1 * s2 * s2 * s2 * s2 * s4 * s4 * s4 +
             B * s0 * s0 * s0 * s0 * s3 * s3 * s4 * s4 * s4 -
             C * s0 * s0 * s0 * s0 * s3 * s3 * s4 * s4 * s4 -
             A * s1 * s1 * s1 * s1 * s3 * s3 * s4 * s4 * s4 +
             C * s1 * s1 * s1 * s1 * s3 * s3 * s4 * s4 * s4 +
             A * s2 * s2 * s2 * s2 * s3 * s3 * s4 * s4 * s4 -
             B * s2 * s2 * s2 * s2 * s3 * s3 * s4 * s4 * s4 -
             B * s0 * s0 * s3 * s3 * s3 * s3 * s4 * s4 * s4 +
             C * s0 * s0 * s3 * s3 * s3 * s3 * s4 * s4 * s4 +
             A * s1 * s1 * s3 * s3 * s3 * s3 * s4 * s4 * s4 -
             C * s1 * s1 * s3 * s3 * s3 * s3 * s4 * s4 * s4 -
             A * s2 * s2 * s3 * s3 * s3 * s3 * s4 * s4 * s4 +
             B * s2 * s2 * s3 * s3 * s3 * s3 * s4 * s4 * s4 -
             C * s0 * s0 * s0 * s1 * s1 * s4 * s4 * s4 * s4 +
             C * s0 * s0 * s1 * s1 * s1 * s4 * s4 * s4 * s4 +
             B * s0 * s0 * s0 * s2 * s2 * s4 * s4 * s4 * s4 -
             A * s1 * s1 * s1 * s2 * s2 * s4 * s4 * s4 * s4 -
             B * s0 * s0 * s2 * s2 * s2 * s4 * s4 * s4 * s4 +
             A * s1 * s1 * s2 * s2 * s2 * s4 * s4 * s4 * s4 -
             B * s0 * s0 * s0 * s3 * s3 * s4 * s4 * s4 * s4 +
             C * s0 * s0 * s0 * s3 * s3 * s4 * s4 * s4 * s4 +
             A * s1 * s1 * s1 * s3 * s3 * s4 * s4 * s4 * s4 -
             C * s1 * s1 * s1 * s3 * s3 * s4 * s4 * s4 * s4 -
             A * s2 * s2 * s2 * s3 * s3 * s4 * s4 * s4 * s4 +
             B * s2 * s2 * s2 * s3 * s3 * s4 * s4 * s4 * s4 +
             B * s0 * s0 * s3 * s3 * s3 * s4 * s4 * s4 * s4 -
             C * s0 * s0 * s3 * s3 * s3 * s4 * s4 * s4 * s4 -
             A * s1 * s1 * s3 * s3 * s3 * s4 * s4 * s4 * s4 +
             C * s1 * s1 * s3 * s3 * s3 * s4 * s4 * s4 * s4 +
             A * s2 * s2 * s3 * s3 * s3 * s4 * s4 * s4 * s4 -
             B * s2 * s2 * s3 * s3 * s3 * s4 * s4 * s4 * s4 -
             D * (s0 - s1) * (s0 - s2) * (s1 - s2) * (s0 - s4) * (s1 - s4) *
                 (s2 - s4) * (s1 * s2 * s4 + s0 * (s2 * s4 + s1 * (s2 + s4)))) /
            ((s0 - s1) * (s0 - s2) * (s1 - s2) * (s0 - s3) * (s1 - s3) *
             (s2 - s3) * (s0 - s4) * (s1 - s4) * (s2 - s4) * (s3 - s4) *
             (s1 * s2 * s3 * s4 +
              s0 * (s2 * s3 * s4 + s1 * (s3 * s4 + s2 * (s3 + s4)))));

  A = f[0].im;
  B = f[1].im;
  C = f[2].im;
  D = f[3].im;
  E = f[4].im;

  c[0].im = (E * s0 * s0 * (s0 - s1) * s1 * s1 * (s0 - s2) * (s1 - s2) * s2 *
                 s2 * (s0 - s3) * (s1 - s3) * (s2 - s3) * s3 * s3 +
             (-D * s0 * s0 * (s0 - s1) * s1 * s1 * (s0 - s2) * (s1 - s2) * s2 *
                  s2 * (s0 - s4) * (s1 - s4) * (s2 - s4) +
              s3 * s3 * (C * s0 * s0 * (s0 - s1) * s1 * s1 * (s0 - s3) *
                             (s1 - s3) * (s0 - s4) * (s1 - s4) +
                         s2 * s2 * (s2 - s3) *
                             (-B * s0 * s0 * (s0 - s2) * (s0 - s3) * (s0 - s4) +
                              A * s1 * s1 * (s1 - s2) * (s1 - s3) * (s1 - s4)) *
                             (s2 - s4)) *
                  (s3 - s4)) *
                 s4 * s4) /
            ((s0 - s1) * (s0 - s2) * (s1 - s2) * (s0 - s3) * (s1 - s3) *
             (s2 - s3) * (s0 - s4) * (s1 - s4) * (s2 - s4) * (s3 - s4) *
             (s1 * s2 * s3 * s4 +
              s0 * (s2 * s3 * s4 + s1 * (s3 * s4 + s2 * (s3 + s4)))));

  c[1].im =
      (-C * s0 * s0 * s0 * s0 * s0 * s1 * s1 * s1 * s1 * s3 * s3 * s3 +
       C * s0 * s0 * s0 * s0 * s1 * s1 * s1 * s1 * s1 * s3 * s3 * s3 +
       B * s0 * s0 * s0 * s0 * s0 * s2 * s2 * s2 * s2 * s3 * s3 * s3 -
       A * s1 * s1 * s1 * s1 * s1 * s2 * s2 * s2 * s2 * s3 * s3 * s3 -
       B * s0 * s0 * s0 * s0 * s2 * s2 * s2 * s2 * s2 * s3 * s3 * s3 +
       A * s1 * s1 * s1 * s1 * s2 * s2 * s2 * s2 * s2 * s3 * s3 * s3 +
       C * s0 * s0 * s0 * s0 * s0 * s1 * s1 * s1 * s3 * s3 * s3 * s3 -
       C * s0 * s0 * s0 * s1 * s1 * s1 * s1 * s1 * s3 * s3 * s3 * s3 -
       B * s0 * s0 * s0 * s0 * s0 * s2 * s2 * s2 * s3 * s3 * s3 * s3 +
       A * s1 * s1 * s1 * s1 * s1 * s2 * s2 * s2 * s3 * s3 * s3 * s3 +
       B * s0 * s0 * s0 * s2 * s2 * s2 * s2 * s2 * s3 * s3 * s3 * s3 -
       A * s1 * s1 * s1 * s2 * s2 * s2 * s2 * s2 * s3 * s3 * s3 * s3 -
       C * s0 * s0 * s0 * s0 * s1 * s1 * s1 * s3 * s3 * s3 * s3 * s3 +
       C * s0 * s0 * s0 * s1 * s1 * s1 * s1 * s3 * s3 * s3 * s3 * s3 +
       B * s0 * s0 * s0 * s0 * s2 * s2 * s2 * s3 * s3 * s3 * s3 * s3 -
       A * s1 * s1 * s1 * s1 * s2 * s2 * s2 * s3 * s3 * s3 * s3 * s3 -
       B * s0 * s0 * s0 * s2 * s2 * s2 * s2 * s3 * s3 * s3 * s3 * s3 +
       A * s1 * s1 * s1 * s2 * s2 * s2 * s2 * s3 * s3 * s3 * s3 * s3 +
       E * (s1 * s1 * s1 * (s1 - s2) * s2 * s2 * s2 * (s1 - s3) * (s2 - s3) *
                s3 * s3 * s3 +
            s0 * s0 * s0 * s0 * s0 *
                (s2 * s2 * s2 * s3 * s3 * s3 * (-s2 + s3) +
                 s1 * s1 * s1 * s1 * (-s2 * s2 * s2 + s3 * s3 * s3) +
                 s1 * s1 * s1 * (s2 * s2 * s2 * s2 - s3 * s3 * s3 * s3)) +
            s0 * s0 * s0 * (s2 * s2 * s2 * s2 * s3 * s3 * s3 * s3 * (-s2 + s3) +
                            s1 * s1 * s1 * s1 * s1 *
                                (-s2 * s2 * s2 * s2 + s3 * s3 * s3 * s3) +
                            s1 * s1 * s1 * s1 * (s2 * s2 * s2 * s2 * s2 -
                                                 s3 * s3 * s3 * s3 * s3)) +
            s0 * s0 * s0 * s0 *
                (s2 * s2 * s2 * s3 * s3 * s3 * (s2 * s2 - s3 * s3) +
                 s1 * s1 * s1 * s1 * s1 * (s2 * s2 * s2 - s3 * s3 * s3) +
                 s1 * s1 * s1 *
                     (-s2 * s2 * s2 * s2 * s2 + s3 * s3 * s3 * s3 * s3))) +
       C * s0 * s0 * s0 * s0 * s0 * s1 * s1 * s1 * s1 * s4 * s4 * s4 -
       C * s0 * s0 * s0 * s0 * s1 * s1 * s1 * s1 * s1 * s4 * s4 * s4 -
       B * s0 * s0 * s0 * s0 * s0 * s2 * s2 * s2 * s2 * s4 * s4 * s4 +
       A * s1 * s1 * s1 * s1 * s1 * s2 * s2 * s2 * s2 * s4 * s4 * s4 +
       B * s0 * s0 * s0 * s0 * s2 * s2 * s2 * s2 * s2 * s4 * s4 * s4 -
       A * s1 * s1 * s1 * s1 * s2 * s2 * s2 * s2 * s2 * s4 * s4 * s4 +
       B * s0 * s0 * s0 * s0 * s0 * s3 * s3 * s3 * s3 * s4 * s4 * s4 -
       C * s0 * s0 * s0 * s0 * s0 * s3 * s3 * s3 * s3 * s4 * s4 * s4 -
       A * s1 * s1 * s1 * s1 * s1 * s3 * s3 * s3 * s3 * s4 * s4 * s4 +
       C * s1 * s1 * s1 * s1 * s1 * s3 * s3 * s3 * s3 * s4 * s4 * s4 +
       A * s2 * s2 * s2 * s2 * s2 * s3 * s3 * s3 * s3 * s4 * s4 * s4 -
       B * s2 * s2 * s2 * s2 * s2 * s3 * s3 * s3 * s3 * s4 * s4 * s4 -
       B * s0 * s0 * s0 * s0 * s3 * s3 * s3 * s3 * s3 * s4 * s4 * s4 +
       C * s0 * s0 * s0 * s0 * s3 * s3 * s3 * s3 * s3 * s4 * s4 * s4 +
       A * s1 * s1 * s1 * s1 * s3 * s3 * s3 * s3 * s3 * s4 * s4 * s4 -
       C * s1 * s1 * s1 * s1 * s3 * s3 * s3 * s3 * s3 * s4 * s4 * s4 -
       A * s2 * s2 * s2 * s2 * s3 * s3 * s3 * s3 * s3 * s4 * s4 * s4 +
       B * s2 * s2 * s2 * s2 * s3 * s3 * s3 * s3 * s3 * s4 * s4 * s4 -
       C * s0 * s0 * s0 * s0 * s0 * s1 * s1 * s1 * s4 * s4 * s4 * s4 +
       C * s0 * s0 * s0 * s1 * s1 * s1 * s1 * s1 * s4 * s4 * s4 * s4 +
       B * s0 * s0 * s0 * s0 * s0 * s2 * s2 * s2 * s4 * s4 * s4 * s4 -
       A * s1 * s1 * s1 * s1 * s1 * s2 * s2 * s2 * s4 * s4 * s4 * s4 -
       B * s0 * s0 * s0 * s2 * s2 * s2 * s2 * s2 * s4 * s4 * s4 * s4 +
       A * s1 * s1 * s1 * s2 * s2 * s2 * s2 * s2 * s4 * s4 * s4 * s4 -
       B * s0 * s0 * s0 * s0 * s0 * s3 * s3 * s3 * s4 * s4 * s4 * s4 +
       C * s0 * s0 * s0 * s0 * s0 * s3 * s3 * s3 * s4 * s4 * s4 * s4 +
       A * s1 * s1 * s1 * s1 * s1 * s3 * s3 * s3 * s4 * s4 * s4 * s4 -
       C * s1 * s1 * s1 * s1 * s1 * s3 * s3 * s3 * s4 * s4 * s4 * s4 -
       A * s2 * s2 * s2 * s2 * s2 * s3 * s3 * s3 * s4 * s4 * s4 * s4 +
       B * s2 * s2 * s2 * s2 * s2 * s3 * s3 * s3 * s4 * s4 * s4 * s4 +
       B * s0 * s0 * s0 * s3 * s3 * s3 * s3 * s3 * s4 * s4 * s4 * s4 -
       C * s0 * s0 * s0 * s3 * s3 * s3 * s3 * s3 * s4 * s4 * s4 * s4 -
       A * s1 * s1 * s1 * s3 * s3 * s3 * s3 * s3 * s4 * s4 * s4 * s4 +
       C * s1 * s1 * s1 * s3 * s3 * s3 * s3 * s3 * s4 * s4 * s4 * s4 +
       A * s2 * s2 * s2 * s3 * s3 * s3 * s3 * s3 * s4 * s4 * s4 * s4 -
       B * s2 * s2 * s2 * s3 * s3 * s3 * s3 * s3 * s4 * s4 * s4 * s4 +
       C * s0 * s0 * s0 * s0 * s1 * s1 * s1 * s4 * s4 * s4 * s4 * s4 -
       C * s0 * s0 * s0 * s1 * s1 * s1 * s1 * s4 * s4 * s4 * s4 * s4 -
       B * s0 * s0 * s0 * s0 * s2 * s2 * s2 * s4 * s4 * s4 * s4 * s4 +
       A * s1 * s1 * s1 * s1 * s2 * s2 * s2 * s4 * s4 * s4 * s4 * s4 +
       B * s0 * s0 * s0 * s2 * s2 * s2 * s2 * s4 * s4 * s4 * s4 * s4 -
       A * s1 * s1 * s1 * s2 * s2 * s2 * s2 * s4 * s4 * s4 * s4 * s4 +
       B * s0 * s0 * s0 * s0 * s3 * s3 * s3 * s4 * s4 * s4 * s4 * s4 -
       C * s0 * s0 * s0 * s0 * s3 * s3 * s3 * s4 * s4 * s4 * s4 * s4 -
       A * s1 * s1 * s1 * s1 * s3 * s3 * s3 * s4 * s4 * s4 * s4 * s4 +
       C * s1 * s1 * s1 * s1 * s3 * s3 * s3 * s4 * s4 * s4 * s4 * s4 +
       A * s2 * s2 * s2 * s2 * s3 * s3 * s3 * s4 * s4 * s4 * s4 * s4 -
       B * s2 * s2 * s2 * s2 * s3 * s3 * s3 * s4 * s4 * s4 * s4 * s4 -
       B * s0 * s0 * s0 * s3 * s3 * s3 * s3 * s4 * s4 * s4 * s4 * s4 +
       C * s0 * s0 * s0 * s3 * s3 * s3 * s3 * s4 * s4 * s4 * s4 * s4 +
       A * s1 * s1 * s1 * s3 * s3 * s3 * s3 * s4 * s4 * s4 * s4 * s4 -
       C * s1 * s1 * s1 * s3 * s3 * s3 * s3 * s4 * s4 * s4 * s4 * s4 -
       A * s2 * s2 * s2 * s3 * s3 * s3 * s3 * s4 * s4 * s4 * s4 * s4 +
       B * s2 * s2 * s2 * s3 * s3 * s3 * s3 * s4 * s4 * s4 * s4 * s4 +
       D * (-s1 * s1 * s1 * (s1 - s2) * s2 * s2 * s2 * (s1 - s4) * (s2 - s4) *
                s4 * s4 * s4 +
            s0 * s0 * s0 * s0 * s0 *
                (s2 * s2 * s2 * (s2 - s4) * s4 * s4 * s4 +
                 s1 * s1 * s1 * s1 * (s2 * s2 * s2 - s4 * s4 * s4) +
                 s1 * s1 * s1 * (-s2 * s2 * s2 * s2 + s4 * s4 * s4 * s4)) +
            s0 * s0 * s0 * s0 *
                (-s2 * s2 * s2 * s2 * s2 * s4 * s4 * s4 +
                 s2 * s2 * s2 * s4 * s4 * s4 * s4 * s4 +
                 s1 * s1 * s1 * s1 * s1 * (-s2 * s2 * s2 + s4 * s4 * s4) +
                 s1 * s1 * s1 *
                     (s2 * s2 * s2 * s2 * s2 - s4 * s4 * s4 * s4 * s4)) +
            s0 * s0 * s0 * (s2 * s2 * s2 * s2 * (s2 - s4) * s4 * s4 * s4 * s4 +
                            s1 * s1 * s1 * s1 * s1 *
                                (s2 * s2 * s2 * s2 - s4 * s4 * s4 * s4) +
                            s1 * s1 * s1 * s1 * (-s2 * s2 * s2 * s2 * s2 +
                                                 s4 * s4 * s4 * s4 * s4)))) /
      ((s0 - s1) * (s0 - s2) * (s1 - s2) * (s0 - s3) * (s1 - s3) * (s2 - s3) *
       (s0 - s4) * (s1 - s4) * (s2 - s4) * (s3 - s4) *
       (s1 * s2 * s3 * s4 +
        s0 * (s2 * s3 * s4 + s1 * (s3 * s4 + s2 * (s3 + s4)))));

  c[2].im =
      (C * s0 * s0 * s0 * s0 * s0 * s1 * s1 * s1 * s1 * s3 * s3 -
       C * s0 * s0 * s0 * s0 * s1 * s1 * s1 * s1 * s1 * s3 * s3 -
       B * s0 * s0 * s0 * s0 * s0 * s2 * s2 * s2 * s2 * s3 * s3 +
       A * s1 * s1 * s1 * s1 * s1 * s2 * s2 * s2 * s2 * s3 * s3 +
       B * s0 * s0 * s0 * s0 * s2 * s2 * s2 * s2 * s2 * s3 * s3 -
       A * s1 * s1 * s1 * s1 * s2 * s2 * s2 * s2 * s2 * s3 * s3 -
       C * s0 * s0 * s0 * s0 * s0 * s1 * s1 * s3 * s3 * s3 * s3 +
       C * s0 * s0 * s1 * s1 * s1 * s1 * s1 * s3 * s3 * s3 * s3 +
       B * s0 * s0 * s0 * s0 * s0 * s2 * s2 * s3 * s3 * s3 * s3 -
       A * s1 * s1 * s1 * s1 * s1 * s2 * s2 * s3 * s3 * s3 * s3 -
       B * s0 * s0 * s2 * s2 * s2 * s2 * s2 * s3 * s3 * s3 * s3 +
       A * s1 * s1 * s2 * s2 * s2 * s2 * s2 * s3 * s3 * s3 * s3 +
       C * s0 * s0 * s0 * s0 * s1 * s1 * s3 * s3 * s3 * s3 * s3 -
       C * s0 * s0 * s1 * s1 * s1 * s1 * s3 * s3 * s3 * s3 * s3 -
       B * s0 * s0 * s0 * s0 * s2 * s2 * s3 * s3 * s3 * s3 * s3 +
       A * s1 * s1 * s1 * s1 * s2 * s2 * s3 * s3 * s3 * s3 * s3 +
       B * s0 * s0 * s2 * s2 * s2 * s2 * s3 * s3 * s3 * s3 * s3 -
       A * s1 * s1 * s2 * s2 * s2 * s2 * s3 * s3 * s3 * s3 * s3 +
       E * (s0 * s0 * s0 * s0 * s0 * (s1 * s1 - s2 * s2) * (s1 * s1 - s3 * s3) *
                (s2 * s2 - s3 * s3) +
            s1 * s1 * s2 * s2 * s3 * s3 *
                (s2 * s2 * s3 * s3 * (-s2 + s3) +
                 s1 * s1 * s1 * (-s2 * s2 + s3 * s3) +
                 s1 * s1 * (s2 * s2 * s2 - s3 * s3 * s3)) +
            s0 * s0 * s0 * s0 *
                (-s2 * s2 * s2 * s2 * s2 * s3 * s3 +
                 s2 * s2 * s3 * s3 * s3 * s3 * s3 +
                 s1 * s1 * s1 * s1 * s1 * (-s2 * s2 + s3 * s3) +
                 s1 * s1 * (s2 * s2 * s2 * s2 * s2 - s3 * s3 * s3 * s3 * s3)) +
            s0 * s0 * (s2 * s2 * s2 * s2 * (s2 - s3) * s3 * s3 * s3 * s3 +
                       s1 * s1 * s1 * s1 * s1 *
                           (s2 * s2 * s2 * s2 - s3 * s3 * s3 * s3) +
                       s1 * s1 * s1 * s1 * (-s2 * s2 * s2 * s2 * s2 +
                                            s3 * s3 * s3 * s3 * s3))) -
       C * s0 * s0 * s0 * s0 * s0 * s1 * s1 * s1 * s1 * s4 * s4 +
       C * s0 * s0 * s0 * s0 * s1 * s1 * s1 * s1 * s1 * s4 * s4 +
       B * s0 * s0 * s0 * s0 * s0 * s2 * s2 * s2 * s2 * s4 * s4 -
       A * s1 * s1 * s1 * s1 * s1 * s2 * s2 * s2 * s2 * s4 * s4 -
       B * s0 * s0 * s0 * s0 * s2 * s2 * s2 * s2 * s2 * s4 * s4 +
       A * s1 * s1 * s1 * s1 * s2 * s2 * s2 * s2 * s2 * s4 * s4 -
       B * s0 * s0 * s0 * s0 * s0 * s3 * s3 * s3 * s3 * s4 * s4 +
       C * s0 * s0 * s0 * s0 * s0 * s3 * s3 * s3 * s3 * s4 * s4 +
       A * s1 * s1 * s1 * s1 * s1 * s3 * s3 * s3 * s3 * s4 * s4 -
       C * s1 * s1 * s1 * s1 * s1 * s3 * s3 * s3 * s3 * s4 * s4 -
       A * s2 * s2 * s2 * s2 * s2 * s3 * s3 * s3 * s3 * s4 * s4 +
       B * s2 * s2 * s2 * s2 * s2 * s3 * s3 * s3 * s3 * s4 * s4 +
       B * s0 * s0 * s0 * s0 * s3 * s3 * s3 * s3 * s3 * s4 * s4 -
       C * s0 * s0 * s0 * s0 * s3 * s3 * s3 * s3 * s3 * s4 * s4 -
       A * s1 * s1 * s1 * s1 * s3 * s3 * s3 * s3 * s3 * s4 * s4 +
       C * s1 * s1 * s1 * s1 * s3 * s3 * s3 * s3 * s3 * s4 * s4 +
       A * s2 * s2 * s2 * s2 * s3 * s3 * s3 * s3 * s3 * s4 * s4 -
       B * s2 * s2 * s2 * s2 * s3 * s3 * s3 * s3 * s3 * s4 * s4 +
       C * s0 * s0 * s0 * s0 * s0 * s1 * s1 * s4 * s4 * s4 * s4 -
       C * s0 * s0 * s1 * s1 * s1 * s1 * s1 * s4 * s4 * s4 * s4 -
       B * s0 * s0 * s0 * s0 * s0 * s2 * s2 * s4 * s4 * s4 * s4 +
       A * s1 * s1 * s1 * s1 * s1 * s2 * s2 * s4 * s4 * s4 * s4 +
       B * s0 * s0 * s2 * s2 * s2 * s2 * s2 * s4 * s4 * s4 * s4 -
       A * s1 * s1 * s2 * s2 * s2 * s2 * s2 * s4 * s4 * s4 * s4 +
       B * s0 * s0 * s0 * s0 * s0 * s3 * s3 * s4 * s4 * s4 * s4 -
       C * s0 * s0 * s0 * s0 * s0 * s3 * s3 * s4 * s4 * s4 * s4 -
       A * s1 * s1 * s1 * s1 * s1 * s3 * s3 * s4 * s4 * s4 * s4 +
       C * s1 * s1 * s1 * s1 * s1 * s3 * s3 * s4 * s4 * s4 * s4 +
       A * s2 * s2 * s2 * s2 * s2 * s3 * s3 * s4 * s4 * s4 * s4 -
       B * s2 * s2 * s2 * s2 * s2 * s3 * s3 * s4 * s4 * s4 * s4 -
       B * s0 * s0 * s3 * s3 * s3 * s3 * s3 * s4 * s4 * s4 * s4 +
       C * s0 * s0 * s3 * s3 * s3 * s3 * s3 * s4 * s4 * s4 * s4 +
       A * s1 * s1 * s3 * s3 * s3 * s3 * s3 * s4 * s4 * s4 * s4 -
       C * s1 * s1 * s3 * s3 * s3 * s3 * s3 * s4 * s4 * s4 * s4 -
       A * s2 * s2 * s3 * s3 * s3 * s3 * s3 * s4 * s4 * s4 * s4 +
       B * s2 * s2 * s3 * s3 * s3 * s3 * s3 * s4 * s4 * s4 * s4 -
       C * s0 * s0 * s0 * s0 * s1 * s1 * s4 * s4 * s4 * s4 * s4 +
       C * s0 * s0 * s1 * s1 * s1 * s1 * s4 * s4 * s4 * s4 * s4 +
       B * s0 * s0 * s0 * s0 * s2 * s2 * s4 * s4 * s4 * s4 * s4 -
       A * s1 * s1 * s1 * s1 * s2 * s2 * s4 * s4 * s4 * s4 * s4 -
       B * s0 * s0 * s2 * s2 * s2 * s2 * s4 * s4 * s4 * s4 * s4 +
       A * s1 * s1 * s2 * s2 * s2 * s2 * s4 * s4 * s4 * s4 * s4 -
       B * s0 * s0 * s0 * s0 * s3 * s3 * s4 * s4 * s4 * s4 * s4 +
       C * s0 * s0 * s0 * s0 * s3 * s3 * s4 * s4 * s4 * s4 * s4 +
       A * s1 * s1 * s1 * s1 * s3 * s3 * s4 * s4 * s4 * s4 * s4 -
       C * s1 * s1 * s1 * s1 * s3 * s3 * s4 * s4 * s4 * s4 * s4 -
       A * s2 * s2 * s2 * s2 * s3 * s3 * s4 * s4 * s4 * s4 * s4 +
       B * s2 * s2 * s2 * s2 * s3 * s3 * s4 * s4 * s4 * s4 * s4 +
       B * s0 * s0 * s3 * s3 * s3 * s3 * s4 * s4 * s4 * s4 * s4 -
       C * s0 * s0 * s3 * s3 * s3 * s3 * s4 * s4 * s4 * s4 * s4 -
       A * s1 * s1 * s3 * s3 * s3 * s3 * s4 * s4 * s4 * s4 * s4 +
       C * s1 * s1 * s3 * s3 * s3 * s3 * s4 * s4 * s4 * s4 * s4 +
       A * s2 * s2 * s3 * s3 * s3 * s3 * s4 * s4 * s4 * s4 * s4 -
       B * s2 * s2 * s3 * s3 * s3 * s3 * s4 * s4 * s4 * s4 * s4 +
       D * (-s0 * s0 * s0 * s0 * s0 * (s1 * s1 - s2 * s2) *
                (s1 * s1 - s4 * s4) * (s2 * s2 - s4 * s4) +
            s1 * s1 * s2 * s2 * s4 * s4 *
                (s2 * s2 * (s2 - s4) * s4 * s4 +
                 s1 * s1 * s1 * (s2 * s2 - s4 * s4) +
                 s1 * s1 * (-s2 * s2 * s2 + s4 * s4 * s4)) +
            s0 * s0 * (s2 * s2 * s2 * s2 * s4 * s4 * s4 * s4 * (-s2 + s4) +
                       s1 * s1 * s1 * s1 * s1 *
                           (-s2 * s2 * s2 * s2 + s4 * s4 * s4 * s4) +
                       s1 * s1 * s1 * s1 *
                           (s2 * s2 * s2 * s2 * s2 - s4 * s4 * s4 * s4 * s4)) +
            s0 * s0 * s0 * s0 *
                (s1 * s1 * s1 * s1 * s1 * (s2 * s2 - s4 * s4) +
                 s2 * s2 * s4 * s4 * (s2 * s2 * s2 - s4 * s4 * s4) +
                 s1 * s1 *
                     (-s2 * s2 * s2 * s2 * s2 + s4 * s4 * s4 * s4 * s4)))) /
      ((s0 - s1) * (s0 - s2) * (s1 - s2) * (s0 - s3) * (s1 - s3) * (s2 - s3) *
       (s0 - s4) * (s1 - s4) * (s2 - s4) * (s3 - s4) *
       (s1 * s2 * s3 * s4 +
        s0 * (s2 * s3 * s4 + s1 * (s3 * s4 + s2 * (s3 + s4)))));

  c[3].im =
      (-C * s0 * s0 * s0 * s0 * s0 * s1 * s1 * s1 * s3 * s3 +
       C * s0 * s0 * s0 * s1 * s1 * s1 * s1 * s1 * s3 * s3 +
       B * s0 * s0 * s0 * s0 * s0 * s2 * s2 * s2 * s3 * s3 -
       A * s1 * s1 * s1 * s1 * s1 * s2 * s2 * s2 * s3 * s3 -
       B * s0 * s0 * s0 * s2 * s2 * s2 * s2 * s2 * s3 * s3 +
       A * s1 * s1 * s1 * s2 * s2 * s2 * s2 * s2 * s3 * s3 +
       C * s0 * s0 * s0 * s0 * s0 * s1 * s1 * s3 * s3 * s3 -
       C * s0 * s0 * s1 * s1 * s1 * s1 * s1 * s3 * s3 * s3 -
       B * s0 * s0 * s0 * s0 * s0 * s2 * s2 * s3 * s3 * s3 +
       A * s1 * s1 * s1 * s1 * s1 * s2 * s2 * s3 * s3 * s3 +
       B * s0 * s0 * s2 * s2 * s2 * s2 * s2 * s3 * s3 * s3 -
       A * s1 * s1 * s2 * s2 * s2 * s2 * s2 * s3 * s3 * s3 -
       C * s0 * s0 * s0 * s1 * s1 * s3 * s3 * s3 * s3 * s3 +
       C * s0 * s0 * s1 * s1 * s1 * s3 * s3 * s3 * s3 * s3 +
       B * s0 * s0 * s0 * s2 * s2 * s3 * s3 * s3 * s3 * s3 -
       A * s1 * s1 * s1 * s2 * s2 * s3 * s3 * s3 * s3 * s3 -
       B * s0 * s0 * s2 * s2 * s2 * s3 * s3 * s3 * s3 * s3 +
       A * s1 * s1 * s2 * s2 * s2 * s3 * s3 * s3 * s3 * s3 +
       E * (s0 * s0 * s0 * s0 * s0 * (s2 * s2 * s3 * s3 * (-s2 + s3) +
                                      s1 * s1 * s1 * (-s2 * s2 + s3 * s3) +
                                      s1 * s1 * (s2 * s2 * s2 - s3 * s3 * s3)) +
            s1 * s1 * s2 * s2 * s3 * s3 *
                (s1 * s1 * s1 * (s2 - s3) + s2 * s3 * (s2 * s2 - s3 * s3) +
                 s1 * (-s2 * s2 * s2 + s3 * s3 * s3)) +
            s0 * s0 * (-s2 * s2 * s2 * s2 * s2 * s3 * s3 * s3 +
                       s2 * s2 * s2 * s3 * s3 * s3 * s3 * s3 +
                       s1 * s1 * s1 * s1 * s1 * (-s2 * s2 * s2 + s3 * s3 * s3) +
                       s1 * s1 * s1 *
                           (s2 * s2 * s2 * s2 * s2 - s3 * s3 * s3 * s3 * s3)) +
            s0 * s0 * s0 * (s1 * s1 * s1 * s1 * s1 * (s2 * s2 - s3 * s3) +
                            s2 * s2 * s3 * s3 * (s2 * s2 * s2 - s3 * s3 * s3) +
                            s1 * s1 * (-s2 * s2 * s2 * s2 * s2 +
                                       s3 * s3 * s3 * s3 * s3))) +
       C * s0 * s0 * s0 * s0 * s0 * s1 * s1 * s1 * s4 * s4 -
       C * s0 * s0 * s0 * s1 * s1 * s1 * s1 * s1 * s4 * s4 -
       B * s0 * s0 * s0 * s0 * s0 * s2 * s2 * s2 * s4 * s4 +
       A * s1 * s1 * s1 * s1 * s1 * s2 * s2 * s2 * s4 * s4 +
       B * s0 * s0 * s0 * s2 * s2 * s2 * s2 * s2 * s4 * s4 -
       A * s1 * s1 * s1 * s2 * s2 * s2 * s2 * s2 * s4 * s4 +
       B * s0 * s0 * s0 * s0 * s0 * s3 * s3 * s3 * s4 * s4 -
       C * s0 * s0 * s0 * s0 * s0 * s3 * s3 * s3 * s4 * s4 -
       A * s1 * s1 * s1 * s1 * s1 * s3 * s3 * s3 * s4 * s4 +
       C * s1 * s1 * s1 * s1 * s1 * s3 * s3 * s3 * s4 * s4 +
       A * s2 * s2 * s2 * s2 * s2 * s3 * s3 * s3 * s4 * s4 -
       B * s2 * s2 * s2 * s2 * s2 * s3 * s3 * s3 * s4 * s4 -
       B * s0 * s0 * s0 * s3 * s3 * s3 * s3 * s3 * s4 * s4 +
       C * s0 * s0 * s0 * s3 * s3 * s3 * s3 * s3 * s4 * s4 +
       A * s1 * s1 * s1 * s3 * s3 * s3 * s3 * s3 * s4 * s4 -
       C * s1 * s1 * s1 * s3 * s3 * s3 * s3 * s3 * s4 * s4 -
       A * s2 * s2 * s2 * s3 * s3 * s3 * s3 * s3 * s4 * s4 +
       B * s2 * s2 * s2 * s3 * s3 * s3 * s3 * s3 * s4 * s4 -
       C * s0 * s0 * s0 * s0 * s0 * s1 * s1 * s4 * s4 * s4 +
       C * s0 * s0 * s1 * s1 * s1 * s1 * s1 * s4 * s4 * s4 +
       B * s0 * s0 * s0 * s0 * s0 * s2 * s2 * s4 * s4 * s4 -
       A * s1 * s1 * s1 * s1 * s1 * s2 * s2 * s4 * s4 * s4 -
       B * s0 * s0 * s2 * s2 * s2 * s2 * s2 * s4 * s4 * s4 +
       A * s1 * s1 * s2 * s2 * s2 * s2 * s2 * s4 * s4 * s4 -
       B * s0 * s0 * s0 * s0 * s0 * s3 * s3 * s4 * s4 * s4 +
       C * s0 * s0 * s0 * s0 * s0 * s3 * s3 * s4 * s4 * s4 +
       A * s1 * s1 * s1 * s1 * s1 * s3 * s3 * s4 * s4 * s4 -
       C * s1 * s1 * s1 * s1 * s1 * s3 * s3 * s4 * s4 * s4 -
       A * s2 * s2 * s2 * s2 * s2 * s3 * s3 * s4 * s4 * s4 +
       B * s2 * s2 * s2 * s2 * s2 * s3 * s3 * s4 * s4 * s4 +
       B * s0 * s0 * s3 * s3 * s3 * s3 * s3 * s4 * s4 * s4 -
       C * s0 * s0 * s3 * s3 * s3 * s3 * s3 * s4 * s4 * s4 -
       A * s1 * s1 * s3 * s3 * s3 * s3 * s3 * s4 * s4 * s4 +
       C * s1 * s1 * s3 * s3 * s3 * s3 * s3 * s4 * s4 * s4 +
       A * s2 * s2 * s3 * s3 * s3 * s3 * s3 * s4 * s4 * s4 -
       B * s2 * s2 * s3 * s3 * s3 * s3 * s3 * s4 * s4 * s4 +
       C * s0 * s0 * s0 * s1 * s1 * s4 * s4 * s4 * s4 * s4 -
       C * s0 * s0 * s1 * s1 * s1 * s4 * s4 * s4 * s4 * s4 -
       B * s0 * s0 * s0 * s2 * s2 * s4 * s4 * s4 * s4 * s4 +
       A * s1 * s1 * s1 * s2 * s2 * s4 * s4 * s4 * s4 * s4 +
       B * s0 * s0 * s2 * s2 * s2 * s4 * s4 * s4 * s4 * s4 -
       A * s1 * s1 * s2 * s2 * s2 * s4 * s4 * s4 * s4 * s4 +
       B * s0 * s0 * s0 * s3 * s3 * s4 * s4 * s4 * s4 * s4 -
       C * s0 * s0 * s0 * s3 * s3 * s4 * s4 * s4 * s4 * s4 -
       A * s1 * s1 * s1 * s3 * s3 * s4 * s4 * s4 * s4 * s4 +
       C * s1 * s1 * s1 * s3 * s3 * s4 * s4 * s4 * s4 * s4 +
       A * s2 * s2 * s2 * s3 * s3 * s4 * s4 * s4 * s4 * s4 -
       B * s2 * s2 * s2 * s3 * s3 * s4 * s4 * s4 * s4 * s4 -
       B * s0 * s0 * s3 * s3 * s3 * s4 * s4 * s4 * s4 * s4 +
       C * s0 * s0 * s3 * s3 * s3 * s4 * s4 * s4 * s4 * s4 +
       A * s1 * s1 * s3 * s3 * s3 * s4 * s4 * s4 * s4 * s4 -
       C * s1 * s1 * s3 * s3 * s3 * s4 * s4 * s4 * s4 * s4 -
       A * s2 * s2 * s3 * s3 * s3 * s4 * s4 * s4 * s4 * s4 +
       B * s2 * s2 * s3 * s3 * s3 * s4 * s4 * s4 * s4 * s4 +
       D * (s1 * s1 * s2 * s2 * s4 * s4 *
                (-s2 * s2 * s2 * s4 + s2 * s4 * s4 * s4 +
                 s1 * s1 * s1 * (-s2 + s4) +
                 s1 * (s2 * s2 * s2 - s4 * s4 * s4)) +
            s0 * s0 * s0 * s0 * s0 *
                (s2 * s2 * (s2 - s4) * s4 * s4 +
                 s1 * s1 * s1 * (s2 * s2 - s4 * s4) +
                 s1 * s1 * (-s2 * s2 * s2 + s4 * s4 * s4)) +
            s0 * s0 * s0 *
                (-s2 * s2 * s2 * s2 * s2 * s4 * s4 +
                 s2 * s2 * s4 * s4 * s4 * s4 * s4 +
                 s1 * s1 * s1 * s1 * s1 * (-s2 * s2 + s4 * s4) +
                 s1 * s1 * (s2 * s2 * s2 * s2 * s2 - s4 * s4 * s4 * s4 * s4)) +
            s0 * s0 * (s2 * s2 * s2 * s4 * s4 * s4 * (s2 * s2 - s4 * s4) +
                       s1 * s1 * s1 * s1 * s1 * (s2 * s2 * s2 - s4 * s4 * s4) +
                       s1 * s1 * s1 * (-s2 * s2 * s2 * s2 * s2 +
                                       s4 * s4 * s4 * s4 * s4)))) /
      ((s0 - s1) * (s0 - s2) * (s1 - s2) * (s0 - s3) * (s1 - s3) * (s2 - s3) *
       (s0 - s4) * (s1 - s4) * (s2 - s4) * (s3 - s4) *
       (s1 * s2 * s3 * s4 +
        s0 * (s2 * s3 * s4 + s1 * (s3 * s4 + s2 * (s3 + s4)))));

  c[4].im = (C * s0 * s0 * s0 * s0 * s1 * s1 * s1 * s3 * s3 -
             C * s0 * s0 * s0 * s1 * s1 * s1 * s1 * s3 * s3 -
             B * s0 * s0 * s0 * s0 * s2 * s2 * s2 * s3 * s3 +
             A * s1 * s1 * s1 * s1 * s2 * s2 * s2 * s3 * s3 +
             B * s0 * s0 * s0 * s2 * s2 * s2 * s2 * s3 * s3 -
             A * s1 * s1 * s1 * s2 * s2 * s2 * s2 * s3 * s3 -
             C * s0 * s0 * s0 * s0 * s1 * s1 * s3 * s3 * s3 +
             C * s0 * s0 * s1 * s1 * s1 * s1 * s3 * s3 * s3 +
             B * s0 * s0 * s0 * s0 * s2 * s2 * s3 * s3 * s3 -
             A * s1 * s1 * s1 * s1 * s2 * s2 * s3 * s3 * s3 -
             B * s0 * s0 * s2 * s2 * s2 * s2 * s3 * s3 * s3 +
             A * s1 * s1 * s2 * s2 * s2 * s2 * s3 * s3 * s3 +
             C * s0 * s0 * s0 * s1 * s1 * s3 * s3 * s3 * s3 -
             C * s0 * s0 * s1 * s1 * s1 * s3 * s3 * s3 * s3 -
             B * s0 * s0 * s0 * s2 * s2 * s3 * s3 * s3 * s3 +
             A * s1 * s1 * s1 * s2 * s2 * s3 * s3 * s3 * s3 +
             B * s0 * s0 * s2 * s2 * s2 * s3 * s3 * s3 * s3 -
             A * s1 * s1 * s2 * s2 * s2 * s3 * s3 * s3 * s3 +
             E * (s0 - s1) * (s0 - s2) * (s1 - s2) * (s0 - s3) * (s1 - s3) *
                 (s2 - s3) * (s1 * s2 * s3 + s0 * (s2 * s3 + s1 * (s2 + s3))) -
             C * s0 * s0 * s0 * s0 * s1 * s1 * s1 * s4 * s4 +
             C * s0 * s0 * s0 * s1 * s1 * s1 * s1 * s4 * s4 +
             B * s0 * s0 * s0 * s0 * s2 * s2 * s2 * s4 * s4 -
             A * s1 * s1 * s1 * s1 * s2 * s2 * s2 * s4 * s4 -
             B * s0 * s0 * s0 * s2 * s2 * s2 * s2 * s4 * s4 +
             A * s1 * s1 * s1 * s2 * s2 * s2 * s2 * s4 * s4 -
             B * s0 * s0 * s0 * s0 * s3 * s3 * s3 * s4 * s4 +
             C * s0 * s0 * s0 * s0 * s3 * s3 * s3 * s4 * s4 +
             A * s1 * s1 * s1 * s1 * s3 * s3 * s3 * s4 * s4 -
             C * s1 * s1 * s1 * s1 * s3 * s3 * s3 * s4 * s4 -
             A * s2 * s2 * s2 * s2 * s3 * s3 * s3 * s4 * s4 +
             B * s2 * s2 * s2 * s2 * s3 * s3 * s3 * s4 * s4 +
             B * s0 * s0 * s0 * s3 * s3 * s3 * s3 * s4 * s4 -
             C * s0 * s0 * s0 * s3 * s3 * s3 * s3 * s4 * s4 -
             A * s1 * s1 * s1 * s3 * s3 * s3 * s3 * s4 * s4 +
             C * s1 * s1 * s1 * s3 * s3 * s3 * s3 * s4 * s4 +
             A * s2 * s2 * s2 * s3 * s3 * s3 * s3 * s4 * s4 -
             B * s2 * s2 * s2 * s3 * s3 * s3 * s3 * s4 * s4 +
             C * s0 * s0 * s0 * s0 * s1 * s1 * s4 * s4 * s4 -
             C * s0 * s0 * s1 * s1 * s1 * s1 * s4 * s4 * s4 -
             B * s0 * s0 * s0 * s0 * s2 * s2 * s4 * s4 * s4 +
             A * s1 * s1 * s1 * s1 * s2 * s2 * s4 * s4 * s4 +
             B * s0 * s0 * s2 * s2 * s2 * s2 * s4 * s4 * s4 -
             A * s1 * s1 * s2 * s2 * s2 * s2 * s4 * s4 * s4 +
             B * s0 * s0 * s0 * s0 * s3 * s3 * s4 * s4 * s4 -
             C * s0 * s0 * s0 * s0 * s3 * s3 * s4 * s4 * s4 -
             A * s1 * s1 * s1 * s1 * s3 * s3 * s4 * s4 * s4 +
             C * s1 * s1 * s1 * s1 * s3 * s3 * s4 * s4 * s4 +
             A * s2 * s2 * s2 * s2 * s3 * s3 * s4 * s4 * s4 -
             B * s2 * s2 * s2 * s2 * s3 * s3 * s4 * s4 * s4 -
             B * s0 * s0 * s3 * s3 * s3 * s3 * s4 * s4 * s4 +
             C * s0 * s0 * s3 * s3 * s3 * s3 * s4 * s4 * s4 +
             A * s1 * s1 * s3 * s3 * s3 * s3 * s4 * s4 * s4 -
             C * s1 * s1 * s3 * s3 * s3 * s3 * s4 * s4 * s4 -
             A * s2 * s2 * s3 * s3 * s3 * s3 * s4 * s4 * s4 +
             B * s2 * s2 * s3 * s3 * s3 * s3 * s4 * s4 * s4 -
             C * s0 * s0 * s0 * s1 * s1 * s4 * s4 * s4 * s4 +
             C * s0 * s0 * s1 * s1 * s1 * s4 * s4 * s4 * s4 +
             B * s0 * s0 * s0 * s2 * s2 * s4 * s4 * s4 * s4 -
             A * s1 * s1 * s1 * s2 * s2 * s4 * s4 * s4 * s4 -
             B * s0 * s0 * s2 * s2 * s2 * s4 * s4 * s4 * s4 +
             A * s1 * s1 * s2 * s2 * s2 * s4 * s4 * s4 * s4 -
             B * s0 * s0 * s0 * s3 * s3 * s4 * s4 * s4 * s4 +
             C * s0 * s0 * s0 * s3 * s3 * s4 * s4 * s4 * s4 +
             A * s1 * s1 * s1 * s3 * s3 * s4 * s4 * s4 * s4 -
             C * s1 * s1 * s1 * s3 * s3 * s4 * s4 * s4 * s4 -
             A * s2 * s2 * s2 * s3 * s3 * s4 * s4 * s4 * s4 +
             B * s2 * s2 * s2 * s3 * s3 * s4 * s4 * s4 * s4 +
             B * s0 * s0 * s3 * s3 * s3 * s4 * s4 * s4 * s4 -
             C * s0 * s0 * s3 * s3 * s3 * s4 * s4 * s4 * s4 -
             A * s1 * s1 * s3 * s3 * s3 * s4 * s4 * s4 * s4 +
             C * s1 * s1 * s3 * s3 * s3 * s4 * s4 * s4 * s4 +
             A * s2 * s2 * s3 * s3 * s3 * s4 * s4 * s4 * s4 -
             B * s2 * s2 * s3 * s3 * s3 * s4 * s4 * s4 * s4 -
             D * (s0 - s1) * (s0 - s2) * (s1 - s2) * (s0 - s4) * (s1 - s4) *
                 (s2 - s4) * (s1 * s2 * s4 + s0 * (s2 * s4 + s1 * (s2 + s4)))) /
            ((s0 - s1) * (s0 - s2) * (s1 - s2) * (s0 - s3) * (s1 - s3) *
             (s2 - s3) * (s0 - s4) * (s1 - s4) * (s2 - s4) * (s3 - s4) *
             (s1 * s2 * s3 * s4 +
              s0 * (s2 * s3 * s4 + s1 * (s3 * s4 + s2 * (s3 + s4)))));
}

#include "PhaseShifts.h"
#include "Basic.h"

delta_params delta_high_energy_tale(delta_params delta) {
  delta_params res;
  double delta_L2, Ddelta_L2; // delta(L2) (d/ds)delta(L2)
  double eps = 1.e-10;

  delta_L2 = gsl_spline_eval(delta.spline, delta.L2, delta.acc);
  Ddelta_L2 = (gsl_spline_eval(delta.spline, delta.L2 + eps, delta.acc) -
               gsl_spline_eval(delta.spline, delta.L2 - eps, delta.acc)) /
              (2. * eps);
  delta.a =
      delta.m * pow(delta.n * M_PI - delta_L2, 2.) / (delta.L2 * Ddelta_L2);
  delta.b = delta.m * (delta.n * M_PI - delta_L2) / (delta.L2 * Ddelta_L2) - 1.;

  res = delta;

  return res;
}

// pipi I=0 S-wave scattering phase shift
double delta0(double s) {
  double v0, v2, sigma0, sigma, k, a0, b0, a2, b2, d0, d2, x0, x1, X0, X1, eps;
  if (s <= delta0_params.s_th) {
    return 0.;
  } else if (s >= delta0_params.L2) {
    return delta0_params.n * M_PI -
           delta0_params.a /
               (delta0_params.b + pow(s / delta0_params.L2, delta0_params.m));
  } else if (pi0pi0 == true) {
    if (s <= 4.) {
      sigma0 = sqrt(1. - 4. * pow(MPION0 / MPION, 2.) / s);
      sigma = sqrt(4. / s - 1.);

      v0 = 0.2212 + 0.2895 * (s - 4.) / 4.;
      v2 = -0.0432 - 0.0803 * (s - 4.) / 4.;

      return atan((v0 + 2. * v2 + 3. * v0 * v2 * sigma) * sigma0 /
                  (3. + 2. * v0 * sigma + v2 * sigma));
    } else {
      sigma = sqrt(1. - 4. / s);
      sigma0 = sqrt(1. - 4. * pow(MPION0 / MPION, 2.) / s);

      d0 = gsl_spline_eval(delta0_params.spline, s, delta0_params.acc);
      d2 = gsl_spline_eval(delta2_params.spline, s, delta2_params.acc);

      x0 = 4.5;
      x1 = 5.5;

      X0 = 10.;
      X1 = 20.;

      if (s < x0) {
        v0 = 0.2212 + 0.2895 * (s - 4.) / 4.;
        v2 = -0.0432 - 0.0803 * (s - 4.) / 4.;

        return atan(
            (2. * v2 * (sigma0 - sigma) +
             v0 * (sigma0 + sigma * (2. + 3. * v2 * v2 * sigma * sigma0))) /
            (3. +
             v2 * sigma *
                 (2. * v0 * (sigma - sigma0) + v2 * (sigma + 2. * sigma0))));
      }

      else if (s < x1) {
        eps = s / (x1 - x0) + x0 / (x0 - x1);

        v0 = eps * tan(d0) / sigma +
             (1. - eps) * (0.2212 + 0.2895 * (s - 4.) / 4.);
        v2 = eps * tan(d2) / sigma +
             (1. - eps) * (-0.0432 - 0.0803 * (s - 4.) / 4.);

        return atan(
            (2. * v2 * (sigma0 - sigma) +
             v0 * (sigma0 + sigma * (2. + 3. * v2 * v2 * sigma * sigma0))) /
            (3. +
             v2 * sigma *
                 (2. * v0 * (sigma - sigma0) + v2 * (sigma + 2. * sigma0))));
      }

      else if (s < X0) {
        v0 = tan(d0) / sigma;
        v2 = tan(d2) / sigma;

        return atan(
            (2. * v2 * (sigma0 - sigma) +
             v0 * (sigma0 + sigma * (2. + 3. * v2 * v2 * sigma * sigma0))) /
            (3. +
             v2 * sigma *
                 (2. * v0 * (sigma - sigma0) + v2 * (sigma + 2. * sigma0))));
      }

      else if (s < X1) {
        eps = s / (X1 - X0) + X0 / (X0 - X1);

        v0 = tan(d0) / sigma;
        v2 = tan(d2) / sigma;

        return eps *
                   gsl_spline_eval(delta0_params.spline, s, delta0_params.acc) +
               (1. - eps) *
                   atan((2. * v2 * (sigma0 - sigma) +
                         v0 * (sigma0 +
                               sigma * (2. + 3. * v2 * v2 * sigma * sigma0))) /
                        (3. +
                         v2 * sigma * (2. * v0 * (sigma - sigma0) +
                                       v2 * (sigma + 2. * sigma0))));
      }

      else {
        return gsl_spline_eval(delta0_params.spline, s, delta0_params.acc);
      }
    }
  } else {
    X0 = 4.1;
    X1 = 4.5;

    if (s <= X0) {
      k = sqrt(s / 4. - 1.);
      return 0.5 * asin(4. * k / sqrt(s) * (0.218581 + 0.287486 * k * k));
    }

    else if (s < X1) {
      eps = s / (X1 - X0) + X0 / (X0 - X1);

      k = sqrt(s / 4. - 1.);

      return eps * gsl_spline_eval(delta0_params.spline, s, delta0_params.acc) +
             (1. - eps) * 0.5 *
                 asin(4. * k / sqrt(s) * (0.218581 + 0.287486 * k * k));
    }

    else {
      return gsl_spline_eval(delta0_params.spline, s, delta0_params.acc);
    }
    // return gsl_spline_eval(delta0_params.spline,s,delta0_params.acc);
  }
}

// pipi I=1 P-wave scattering phase shift
double delta1(double s) {
  if (s <= delta1_params.s_th) {
    return 0.;
  } else if (s >= delta1_params.L2) {
    return delta1_params.n * M_PI -
           delta1_params.a /
               (delta1_params.b + pow(s / delta1_params.L2, delta1_params.m));
  } else {
    return gsl_spline_eval(delta1_params.spline, s, delta1_params.acc);
  }
}

// Gibt Phase der P-Welle (pipi-Streuung) in Abhängigkeit der Mandelstamvariable
// s wieder. Parametrisierung aus: https://arxiv.org/pdf/1102.2183.pdf
/*double delta1(double s) {
    s *= MPION*MPION;
    double Mk=493.7;
    double c=0;
    double f1,f2;
    double Mrho=773.6, Mp=139.6, B0=1.055, B1=0.15, s00=1050*1050; //Konstanten
für f1
    double k=sqrt(s/4.0-Mp*Mp);
    double aa=sqrt(s)/(2*k*k*k)*(Mrho*Mrho-s);
    double bb=2*Mp*Mp*Mp/(Mrho*Mrho*sqrt(s));
    double cc=B0;
    double dd=B1*(sqrt(s)-sqrt(s00-s))/(sqrt(s)+sqrt(s00-s));
    double a0=2.675093, a1=1.39, a2=-1.70; //Lambda im paper; Konstanten für f2
    double a=-388199.429343, b=288200.109581;
    double s0=1742400.000000;
    if (s<=4*Mp*Mp) return 0;
    if (s<=4*Mk*Mk) {
        f1=atan(1.0/(aa*(bb+cc+dd)));
        if (f1<0) c=M_PI;
        else c=0;
        return f1+c;
    }
    else if (s<=s0) {
        f2=a0+a1*(sqrt(s)/(2.0*Mk)-1)+a2*(sqrt(s)/(2.0*Mk)-1)*(sqrt(s)/(2.0*Mk)-1);
        return f2;
    }
    else return M_PI+a/(b+s);
}*/

// pipi I=2 S-wave scattering phase shift
double delta2(double s) {
  if (s <= delta2_params.s_th) {
    return 0.;
  }
  // else if (s>=1000.) {
  //    return 0.;
  //}
  else if (s >= delta2_params.L2) {
    return delta2_params.n * M_PI -
           delta2_params.a /
               (delta2_params.b + pow(s / delta2_params.L2, delta2_params.m));
  } else {
    return gsl_spline_eval(delta2_params.spline, s, delta2_params.acc);
  }
}

// etapi I=1 S-wave scattering phase shift
double delta_etapi(double s) {
  if (s <= delta_etapi_params.s_th) {
    return 0.;
  } else if (s >= delta_etapi_params.L2) {
    return delta_etapi_params.n * M_PI -
           delta_etapi_params.a /
               (delta_etapi_params.b +
                pow(s / delta_etapi_params.L2, delta_etapi_params.m));
  } else if (pi0pi0 == true) {
    if (s < 4. * pow(MKAON / MPION, 2.)) {
      double m, n, sp;

      m = (4. * pow(MKAON / MPION, 2.) - pow(1. + META / MPION, 2.)) /
          (4. * pow(MKAON / MPION, 2.) -
           pow(MPION0 / MPION + META / MPION, 2.));
      n = 4. * pow(MKAON / MPION, 2.) * (1. - m);

      return gsl_spline_eval(delta_etapi_params.spline, m * s + n,
                             delta_etapi_params.acc);
    } else {
      return gsl_spline_eval(delta_etapi_params.spline, s,
                             delta_etapi_params.acc);
    }
  } else {
    return gsl_spline_eval(delta_etapi_params.spline, s,
                           delta_etapi_params.acc);
  }
}

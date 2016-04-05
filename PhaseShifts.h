#ifndef _PHASESHIFTS_H_
#define _PHASESHIFTS_H_

gsl_spline *delta0_spline; //pipi I=0 S-wave scattering phase shift
gsl_spline *delta1_spline; //pipi I=1 P-wave scattering phase shift
gsl_spline *delta2_spline; //pipi I=2 S-wave scattering phase shift
gsl_spline *delta_etapi_spline; //etapi I=1 S-wave scattering phase shift

gsl_interp_accel *acc_delta0; //gsl spline accelerator for delta0
gsl_interp_accel *acc_delta1; //gsl spline accelerator for delta1
gsl_interp_accel *acc_delta2; //gsl spline accelerator for delta2
gsl_interp_accel *acc_delta_etapi; //gsl spline accelerator for delta_etapi

double delta0_const; //value of delta0 at cutoff L2
double delta1_const; //value of delta1 at cutoff L2
double delta2_const; //value of delta2 at cutoff L2
double delta_etapi_const; //value of delta_etapi at cutoff L2

double L2_delta0; //cutoff L2 for delta0
double L2_delta1; //cutoff L2 for delta1
double L2_delta2; //cutoff L2 for delta2
double L2_delta_etapi; //cutoff L2 for delta_etapi

//pipi I=0 S-wave scattering phase shift
double delta0(double s);

//pipi I=1 P-wave scattering phase shift
double delta1(double s);

//pipi I=2 S-wave scattering phase shift
double delta2(double s);

//etapi I=1 S-wave scattering phase shift
double delta_etapi(double s);

#endif

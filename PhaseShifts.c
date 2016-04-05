#include "Basic.h"
#include "PhaseShifts.h"

//pipi I=0 S-wave scattering phase shift
double delta0(double s) {
    if (s<=s_pipi) {
        return 0.;
    }
    else if (s>=L2_delta0) {
        return delta0_const;
    }
    else {
        return gsl_spline_eval(delta0_spline,s,acc_delta0);
    }
}

//pipi I=1 P-wave scattering phase shift
double delta1(double s) {
    if (s<=s_pipi) {
        return 0.;
    }
    else if (s>=L2_delta1) {
        return delta1_const;
    }
    else {
        return gsl_spline_eval(delta1_spline,s,acc_delta1);
    }
}

//pipi I=2 S-wave scattering phase shift
double delta2(double s) {
    if (s<=s_pipi) {
        return 0.;
    }
    else if (s>=L2_delta2) {
        return delta2_const;
    }
    else {
        return gsl_spline_eval(delta2_spline,s,acc_delta2);
    }
}

//etapi I=1 S-wave scattering phase shift
double delta_etapi(double s) {
    if (s<=s_etapi) {
        return 0.;
    }
    else if (s>=L2_delta_etapi) {
        return delta_etapi_const;
    }
    else {
        return gsl_spline_eval(delta_etapi_spline,s,acc_delta_etapi);
    }
}

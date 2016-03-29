#include "Basic.h"
#include "PhaseShifts.h"

//paramterisation region end of the phase shifts in MPION^2
#define S_END 90.0

//parameters for schenk phase shifts
//delta00
#define SCHENK_S00 36.77
#define SCHENK_A00 0.220
#define SCHENK_B00 0.268
#define SCHENK_C00 -0.139e-01
#define SCHENK_D00 -0.139e-02
#define SCHENK_CONST00  1.5*M_PI
#define SCHENK_KS_A00 4.5963643469 //delta00_schenk(S_END)
#define SCHENK_KS_B00 0.0024247271 //(d/ds)delta00_schenk(S_END))
const double SCHENK_KS_C00 = 3.0*(SCHENK_CONST00-SCHENK_KS_A00-2.0/3.0*SCHENK_KS_B00*(LAMBDA_SQUARED-S_END))/((LAMBDA_SQUARED-S_END)*(LAMBDA_SQUARED-S_END));
const double SCHENK_KS_D00 = -(SCHENK_KS_B00+6.0*(SCHENK_CONST00-SCHENK_KS_A00-2.0/3.0*SCHENK_KS_B00*(LAMBDA_SQUARED-S_END))/(LAMBDA_SQUARED-S_END))/(3.0*(LAMBDA_SQUARED-S_END)*(LAMBDA_SQUARED-S_END));

//delta11
#define SCHENK_S11 30.72
#define SCHENK_A11 0.379e-01
#define SCHENK_B11 0.140e-04
#define SCHENK_C11 -0.673e-04
#define SCHENK_D11 0.163e-07
#define SCHENK_CONST11 M_PI
#define SCHENK_KS_A11 3.0729862911 //delta11_schenk(S_END)
#define SCHENK_KS_B11 0.0070743411 //(d/ds)delta11_schenk(S_END))
const double SCHENK_KS_C11 = 3.0*(SCHENK_CONST11-SCHENK_KS_A11-2.0/3.0*SCHENK_KS_B11*(LAMBDA_SQUARED-S_END))/((LAMBDA_SQUARED-S_END)*(LAMBDA_SQUARED-S_END));
const double SCHENK_KS_D11 = -(SCHENK_KS_B11+6.0*(SCHENK_CONST11-SCHENK_KS_A11-2.0/3.0*SCHENK_KS_B11*(LAMBDA_SQUARED-S_END))/(LAMBDA_SQUARED-S_END))/(3.0*(LAMBDA_SQUARED-S_END)*(LAMBDA_SQUARED-S_END));

//delta20
#define SCHENK_S20 -21.62
#define SCHENK_A20 -0.444e-01
#define SCHENK_B20 -0.857e-01
#define SCHENK_C20 -0.221e-02
#define SCHENK_D20 -0.129e-03
#define SCHENK_CONST20 -0.5*M_PI
#define SCHENK_KS_A20 -1.5684955005 //delta20_schenk(S_END)
#define SCHENK_KS_B20 -0.0000821565 //(d/ds)delta20_schenk(S_END))
const double SCHENK_KS_C20 = 3.0*(SCHENK_CONST20-SCHENK_KS_A20-2.0/3.0*SCHENK_KS_B20*(LAMBDA_SQUARED-S_END))/((LAMBDA_SQUARED-S_END)*(LAMBDA_SQUARED-S_END));
const double SCHENK_KS_D20 = -(SCHENK_KS_B20+6.0*(SCHENK_CONST20-SCHENK_KS_A20-2.0/3.0*SCHENK_KS_B20*(LAMBDA_SQUARED-S_END))/(LAMBDA_SQUARED-S_END))/(3.0*(LAMBDA_SQUARED-S_END)*(LAMBDA_SQUARED-S_END));

//parameters for madrid phase shifts
//delta00
#define MADRID_00_B0 7.26
#define MADRID_00_B1 -25.3
#define MADRID_00_B2 -33.1
#define MADRID_00_B3 -26.6
#define MADRID_00_D0 227.1
#define MADRID_00_C0 660.0
#define MADRID_00_B 94.0
#define MADRID_00_C 40.4
#define MADRID_00_D -86.9
#define MADRID_00_deltaM 92.2069882491 //delta00_madrid(sM)
#define MADRID_00_deltaMp 1.8022295047 //(d/ds)delta00_madrid(sM)
const double MADRID_00_SM = (850.0*850.0)/(MPION*MPION);
#define MADRID_CONST00  2.0*M_PI
#define MADRID_KS_A00 4.9962715201 //delta00_madrid(S_END)
#define MADRID_KS_B00 0.0301803027 //(d/ds)delta00_madrid(S_END)
const double MADRID_KS_C00 = 3.0*(MADRID_CONST00-MADRID_KS_A00-2.0/3.0*MADRID_KS_B00*(LAMBDA_SQUARED-S_END))/((LAMBDA_SQUARED-S_END)*(LAMBDA_SQUARED-S_END));
const double MADRID_KS_D00 = -(MADRID_KS_B00+6.0*(MADRID_CONST00-MADRID_KS_A00-2.0/3.0*MADRID_KS_B00*(LAMBDA_SQUARED-S_END))/(LAMBDA_SQUARED-S_END))/(3.0*(LAMBDA_SQUARED-S_END)*(LAMBDA_SQUARED-S_END));

//delta11
#define MADRID_11_B0 1.055
#define MADRID_11_B1 0.15
#define MADRID_11_LAMBDA0 2.6809640124 //delta11_madrid(4.0*MKAON^2)
#define MADRID_11_LAMBDA1 1.57
#define MADRID_11_LAMBDA2 -1.96
#define MADRID_CONST11  M_PI
#define MADRID_KS_A11 2.9879168835 //delta11_madrid(S_END)
#define MADRID_KS_B11 0.0017963409 //(d/ds)delta11_madrid(S_END))
const double MADRID_11_S0 = (1050.0*1050.0)/(MPION*MPION);
const double MADRID_KS_C11 = 3.0*(MADRID_CONST11-MADRID_KS_A11-2.0/3.0*MADRID_KS_B11*(LAMBDA_SQUARED-S_END))/((LAMBDA_SQUARED-S_END)*(LAMBDA_SQUARED-S_END));
const double MADRID_KS_D11 = -(MADRID_KS_B11+6.0*(MADRID_CONST11-MADRID_KS_A11-2.0/3.0*MADRID_KS_B11*(LAMBDA_SQUARED-S_END))/(LAMBDA_SQUARED-S_END))/(3.0*(LAMBDA_SQUARED-S_END)*(LAMBDA_SQUARED-S_END));

//delta20
#define MADRID_20_B0 -80.4
#define MADRID_20_B1 -73.6
#define MADRID_20_WLM 0.1592689620 //W-Function(sM,sl)
#define MADRID_20_WHM -0.1446529057 //W-Function(sM,sh)
#define MADRID_20_BH2 112.0
const double MADRID_20_BH0 = MADRID_20_B0+MADRID_20_B1*MADRID_20_WLM;
const double MADRID_20_BH1 = MADRID_20_B1*1.8532931282;
const double MADRID_20_SM = (850.0*850.0)/(MPION*MPION);
const double MADRID_20_SL = (1050.0*1050.0)/(MPION*MPION);
const double MADRID_20_SH = (1420.0*1420.0)/(MPION*MPION);

#define MADRID_CONST20  0.0
#define MADRID_KS_A20 -0.5719473698 //delta20_madrid(S_END)
#define MADRID_KS_B20 -0.0051114668 //(d/ds)delta20_madrid(S_END))
const double MADRID_KS_C20 = 3.0*(MADRID_CONST20-MADRID_KS_A20-2.0/3.0*MADRID_KS_B20*(LAMBDA_SQUARED-S_END))/((LAMBDA_SQUARED-S_END)*(LAMBDA_SQUARED-S_END));
const double MADRID_KS_D20 = -(MADRID_KS_B20+6.0*(MADRID_CONST20-MADRID_KS_A20-2.0/3.0*MADRID_KS_B20*(LAMBDA_SQUARED-S_END))/(LAMBDA_SQUARED-S_END))/(3.0*(LAMBDA_SQUARED-S_END)*(LAMBDA_SQUARED-S_END));

//Pion-Pion Scattering Phase Shifts deltaIJ (I: Isospin, J: Angular Momentum) given for s in MPION^2

double kubic_spline(double s, double sCut, double sConst, double A, double B, double fconst) {
    double C,D,ds;
    
    ds = sConst-sCut;
    
    C = 3.0*(fconst-A-2.0/3.0*B*ds)/pow(ds,2.0);
    D = -(B+2.0*C*ds)/(3.0*pow(ds,2.0));
    
    ds = s-sCut;
    
    return A+B*ds+C*pow(ds,2.0)+D*pow(ds,3.0);
}

//Schenk Parametrisation

double delta00_schenk(double s) {
    double q2,result;
    
    if (s<=4.0) {
        result = 0.0;
    }
    
    else if (s>4.0 && s<SCHENK_S00) {
        q2 = Q_SQUARE(s);
        result = atan(sqrt(SIGMA_SQUARE(s))*(4.0-SCHENK_S00)/(s-SCHENK_S00)*(SCHENK_A00+SCHENK_B00*q2+SCHENK_C00*pow(q2,2.0)+SCHENK_D00*pow(q2,3.0)));
    }
    
    else if (s==SCHENK_S00) {
        result = 0.5*M_PI;
    }
    
    else if (s>SCHENK_S00 && s<=S_END) {
        q2 = Q_SQUARE(s);
        result = atan(sqrt(SIGMA_SQUARE(s))*(4.0-SCHENK_S00)/(s-SCHENK_S00)*(SCHENK_A00+SCHENK_B00*q2+SCHENK_C00*pow(q2,2.0)+SCHENK_D00*pow(q2,3.0)))+M_PI;
    }
    
    else if (s>S_END && s<=LAMBDA_SQUARED) {
        q2 = (s-S_END);
        
        result = SCHENK_KS_A00+SCHENK_KS_B00*q2+SCHENK_KS_C00*pow(q2,2.0)+SCHENK_KS_D00*pow(q2,3.0);
    }
    
    else if (s>LAMBDA_SQUARED) {
        result = SCHENK_CONST00;
    }
    
    return result;
}

double delta11_schenk(double s) {
    double q2,result;
    
    if (s<=4.0) {
        result = 0.0;
    }
    
    else if (s>4.0 && s<SCHENK_S11) {
        q2 = Q_SQUARE(s);
        result = atan(sqrt(SIGMA_SQUARE(s))*q2*(4.0-SCHENK_S11)/(s-SCHENK_S11)*(SCHENK_A11+SCHENK_B11*q2+SCHENK_C11*q2*q2+SCHENK_D11*pow(q2,3.0)));
    }
    
    else if (s==SCHENK_S11) {
        result = M_PI/2.0;
    }
    
    else if (s>SCHENK_S11 && s<=S_END) {
        q2 = Q_SQUARE(s);
        result = atan(sqrt(SIGMA_SQUARE(s))*q2*(4.0-SCHENK_S11)/(s-SCHENK_S11)*(SCHENK_A11+SCHENK_B11*q2+SCHENK_C11*q2*q2+SCHENK_D11*pow(q2,3.0)))+M_PI;
    }
    
    else if (s>S_END && s<=LAMBDA_SQUARED) {
        q2 = (s-S_END);
        
        result = SCHENK_KS_A11+SCHENK_KS_B11*q2+SCHENK_KS_C11*pow(q2,2.0)+SCHENK_KS_D11*pow(q2,3.0);
    }
    
    else if (s>LAMBDA_SQUARED) {
        result = SCHENK_CONST11;
    }
    
    return result;
}

double delta20_schenk(double s) {
    double q2,result;
    
    if (s<=4.0) {
        result = 0.0;
    }
    
    if (s>4.0 && s<=S_END) {
        q2 = Q_SQUARE(s);
        result = atan(sqrt(SIGMA_SQUARE(s))*q2*q2*(4.0-SCHENK_S20)/(s-SCHENK_S20)*(SCHENK_A20+SCHENK_B20*q2+SCHENK_C20*q2*q2+SCHENK_D20*pow(q2,3.0)));
    }
    
    else if (s>S_END && s<=LAMBDA_SQUARED) {
        q2 = (s-S_END);
        
        result = SCHENK_KS_A20+SCHENK_KS_B20*q2+SCHENK_KS_C20*pow(q2,2.0)+SCHENK_KS_D20*pow(q2,3.0);
    }
    
    else if (s>LAMBDA_SQUARED) {
        result = SCHENK_CONST20;
    }
    
    return result;
}

//Madrid Parametrisation

double delta00_madrid(double s) {
    double q,k,result;
    
    //p = (B-2.0*C+D*pow(MKAON/META,2.0))*4.0*pow(MKAON,2.0)/C;
    //q = (d0-B+C-D-360.0)*16.0*pow(MKAON,4.0)/C;
    //s2pi = (-p/2.0+sqrt(pow(p/2.0,2.0)-q))/pow(MPION,2.0);
    
    if (s<=4.0) {
        result = 0.0;
    }
    
    else if (s>4.0 && s<=MADRID_00_SM) {
        k = sqrt(Q_SQUARE(s));
        q = W_FUNCTION(s,4.0*MKAON_SQUARED);
        
        result = atan(2.0*k*(s-0.5)/(sqrt(s)*(1.0/sqrt(s)+MADRID_00_B0+MADRID_00_B1*q+MADRID_00_B2*pow(q,2.0)+MADRID_00_B3*pow(q,3.0))));
        if (result<0.0) {
            result += M_PI;
        }
    }
    
    else if (s>MADRID_00_SM && s<4.0*MKAON_SQUARED) {
        k = sqrt(-Q1_SQUARE(s));
        q = sqrt(-Q1_SQUARE(MADRID_00_SM));
        
        result = (MADRID_00_D0*pow(1.0-k/q,2.0)+MADRID_00_deltaM*k/q*(2.0-k/q)+k*(q-k)*(8.0*MADRID_00_deltaMp+MADRID_00_C0*(q-k)/pow(MKAON_SQUARED,1.5)))*M_PI/180.0;
    }
    
    else if (s>=4.0*MKAON_SQUARED && s<=4.0*META_SQUARED) {
        k = Q1_SQUARE(s)/MKAON_SQUARED;
        
        result = (MADRID_00_D0+MADRID_00_B*k+MADRID_00_C*pow(k,2.0))*M_PI/180.0;
    }
    
    else if (s>4.0*META_SQUARED && s<=S_END) {
        k = Q1_SQUARE(s)/MKAON_SQUARED;
        q = Q2_SQUARE(s)/META_SQUARED;
        
        result = (MADRID_00_D0+MADRID_00_B*k+MADRID_00_C*pow(k,2.0)+MADRID_00_D*q)*M_PI/180.0;
    }
    
    else if (s>S_END && s<=LAMBDA_SQUARED) {
        k = (s-S_END);
        
        result = MADRID_KS_A00+MADRID_KS_B00*k+MADRID_KS_C00*pow(k,2.0)+MADRID_KS_D00*pow(k,3.0);
    }
    
    else if (s>LAMBDA_SQUARED) {
        result = MADRID_CONST00;
    }
    
    return result;
}

double delta11_madrid(double s) {
    double k,w,result;
    
    if (s<=4.0) {
        result = 0.0;
    }
    
    else if (s>4.0 && s<=4.0*MKAON_SQUARED) {
        k = Q_SQUARE(s);
        k *= sqrt(k);
        w = W_FUNCTION(s,MADRID_11_S0);
        
        result = atan(2.0*k/(sqrt(s)*(MRHO_SQUARED-s)*(2.0/(MRHO_SQUARED*sqrt(s))+MADRID_11_B0+MADRID_11_B1*w)));
        
        if (result<0.0) {
            result += M_PI;
        }
    }
    
    else if (s>4.0*MKAON_SQUARED && s<=S_END) {
        result = MADRID_11_LAMBDA0+MADRID_11_LAMBDA1*(0.5*sqrt(s/MKAON_SQUARED)-1.0)+MADRID_11_LAMBDA2*pow(0.5*sqrt(s/MKAON_SQUARED)-1.0,2.0);
    }
    
    else if (s>S_END && s<=LAMBDA_SQUARED) {
        k = (s-S_END);
        
        result = MADRID_KS_A11+MADRID_KS_B11*k+MADRID_KS_C11*pow(k,2.0)+MADRID_KS_D11*pow(k,3.0);
    }
    
    else if (s>LAMBDA_SQUARED) {
        result = MADRID_CONST11;
    }
    
    return result;
}

double delta20_madrid(double s) {
    double k,w,result;
    
    if (s<=4.0) {
        result = 0.0;
    }
    
    else if (s>4.0 && s<=MADRID_20_SM) {
        k = sqrt(Q_SQUARE(s));
        w = W_FUNCTION(s,MADRID_20_SL);
        
        result = atan(2.0*k*(s-2.0)/(sqrt(s)*(MADRID_20_B0+MADRID_20_B1*w)));
    }
    
    else if (s>MADRID_20_SM && s<=S_END) {
        w = W_FUNCTION(s,MADRID_20_SH);
        k = sqrt(Q_SQUARE(s));
        
        result = atan(2.0*k*(s-2.0)/(sqrt(s)*(MADRID_20_BH0+MADRID_20_BH1*(w-MADRID_20_WHM)+MADRID_20_BH2*pow(w-MADRID_20_WHM,2.0))));
    }
    
    else if (s>S_END && s<=LAMBDA_SQUARED) {
        k = (s-S_END);
        
        result = MADRID_KS_A20+MADRID_KS_B20*k+MADRID_KS_C20*pow(k,2.0)+MADRID_KS_D20*pow(k,3.0);
    }
    
    else if (s>LAMBDA_SQUARED) {
        result = MADRID_CONST20;
    }
    
    return result;
}

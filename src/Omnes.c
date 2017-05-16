#include "Basic.h"
#include "AngularAverages.h"
#include "Omnes.h"

double dOmnes0(double (*delta)(double), delta_params Delta) {
    int M = 1500;
    double eps = 1.e-7,res,tmp,err;
    gsl_integration_workspace *w = gsl_integration_workspace_alloc(M);
    gsl_function F;
    
    double f(double z, void * params) {
        double fct = delta(z)/(z*z);
        return fct;
    }
    
    double high_energy_reminder(double z, void * params) {
        double fct = Delta.a/(Delta.L2*(Delta.b+pow(1.-z,-Delta.m)));
        return fct;
    }
    
    F.function = &f;
    F.params = NULL;
    
    gsl_integration_qags(&F,Delta.s_th,Delta.L2,eps,eps,M,w,&res,&err);
    
    //printf("threshold: %.3e\n",res);
    
    F.function = &high_energy_reminder;
    F.params = NULL;
    
    gsl_integration_qags(&F,0.,1.,eps,eps,M,w,&tmp,&err);
    res += Delta.n*M_PI/Delta.L2-tmp;
    
    //printf("infty: %.3e\n",Delta.n*M_PI/Delta.L2-tmp);
    
    return res/M_PI;
}

double d2Omnes0(double (*delta)(double), delta_params Delta, double dO0) {
    int M = 1500;
    double eps = 1.e-7,res,tmp,err;
    gsl_integration_workspace *w = gsl_integration_workspace_alloc(M);
    gsl_function F;
    
    double f(double z, void * params) {
        double fct = delta(z)/(z*z*z);
        return fct;
    }
    
    double high_energy_reminder(double z, void * params) {
        double fct = Delta.a*(1.-z)/(Delta.L2*Delta.L2*(Delta.b+pow(1.-z,-Delta.m)));
        return fct;
    }
    
    F.function = &f;
    F.params = NULL;
    
    gsl_integration_qags(&F,Delta.s_th,Delta.L2,eps,eps,M,w,&res,&err);
    
    F.function = &high_energy_reminder;
    F.params = NULL;
    
    gsl_integration_qags(&F,0.,1.,eps,eps,M,w,&tmp,&err);
    res += 0.5*Delta.n*M_PI/pow(Delta.L2,2.)-tmp;
    
    return 2.*res/M_PI+dO0*dO0;
}

complex Omnes_function(double (*delta)(double), delta_params Delta, complex s) {
    int M = 1500;
    complex temp;
    double A,B,C,D,r,phi,result,result1,result2,error,eps,eps1,eps2,s_temp;
    double pts[3];
    
    gsl_integration_workspace *w = gsl_integration_workspace_alloc(M);
    gsl_function F;
    
    eps = 1.0e-7;
    eps1 = 1.0e-07;
    eps2 = 1.0e-07;
    
    pts[0] = Delta.s_th;
    pts[2] = Delta.L2;
    
    if (fabs(s.im)<eps) {
        double f1 (double z, void * params) {
            double x = *(double *) params;
            double fct = delta(z)/(z*(z-x));
            return fct;
        }
        
        double f2 (double z, void * params) {
            double x = *(double *) params;
            double fct = (delta(z)-delta(x))/(z*(z-x));
            return fct;
        }
        
        double high_energy_reminder(double z, void * params) {
            double x = *(double *)params;
            double fct = Delta.a*pow(1.-z,Delta.m)/((Delta.b*pow(1.-z,Delta.m)+1.)*(Delta.L2-x*(1.-z)));
            return fct;
        }
        
        double high_energy_reminder_cv(double z, void * params) {
            double y = *(double *)params;
            double fct = Delta.a*pow(1.-z,Delta.m)*(1.-y)/((Delta.b*pow(1.-z,Delta.m)+1.)*Delta.L2);
            return fct;
        }
        
        if (s.re<Delta.s_th-eps) {
            F.function = &f1;
            F.params = &s.re;
            
            gsl_integration_qags(&F,Delta.s_th,Delta.L2,eps1,eps2,M,w,&result,&error);
            
            F.function = &high_energy_reminder;
            F.params = &s.re;
            
            gsl_integration_qags(&F,0.,1.,eps1,eps2,M,w,&result1,&error);
            
            result = s.re/M_PI*(result-result1)+Delta.n*log(Delta.L2/(Delta.L2-s.re));
        }
        
        else if (fabs(s.re-Delta.s_th)<eps) {
            F.function = &f1;
            s_temp = Delta.s_th-eps;
            F.params = &s_temp;
            
            gsl_integration_qags(&F,Delta.s_th,Delta.L2,eps1,eps2,M,w,&result,&error);
            
            F.function = &high_energy_reminder;
            F.params = &s_temp;
            
            gsl_integration_qags(&F,0.,1.,eps1,eps2,M,w,&result1,&error);
            
            result = s.re/M_PI*(result-result1)+Delta.n*log(Delta.L2/(Delta.L2-s.re));
            
            F.function = &f2;
            s_temp = Delta.s_th+eps;
            F.params = &s_temp;
            pts[1] = s_temp;
            
            gsl_integration_qagp(&F,pts,3,eps1,eps2,M,w,&result1,&error);
            
            F.function = &high_energy_reminder;
            F.params = &s_temp;
            
            gsl_integration_qags(&F,0.,1.,eps1,eps2,M,w,&result2,&error);
            
            result1 = ((result1-result2)*s.re+(log(Delta.s_th/(s_temp-Delta.s_th))+log((Delta.L2-s.re)/Delta.L2))*delta(s.re))/M_PI+Delta.n*log(Delta.L2/(Delta.L2-s.re));
            
            result = 0.5*(result+result1);
        }
        
        else if (s.re>Delta.s_th+eps && s.re<Delta.L2-eps) {
            F.function = &f2;
            F.params = &s.re;
            pts[1] = s.re;
            
            gsl_integration_qagp(&F,pts,3,eps1,eps2,M,w,&result,&error);
            
            F.function = &high_energy_reminder;
            F.params = &s.re;
            
            gsl_integration_qags(&F,0.,1.,eps1,eps2,M,w,&result1,&error);
            
            result = ((result-result1)*s.re+(log(Delta.s_th/(s.re-Delta.s_th))+log((Delta.L2-s.re)/Delta.L2))*delta(s.re))/M_PI+Delta.n*log(Delta.L2/(Delta.L2-s.re));
        }
        
        else if (fabs(s.re-Delta.L2)<eps) {
            F.function = &f1;
            s_temp = Delta.L2+eps;
            F.params = &s_temp;
            pts[1] = s_temp;
            
            gsl_integration_qagp(&F,pts,3,eps1,eps2,M,w,&result,&error);
            
            F.function = &high_energy_reminder;
            F.params = &s_temp;
            
            gsl_integration_qags(&F,0.,1.,eps1,eps2,M,w,&result1,&error);
            
            result = s.re/M_PI*(result-result1)+Delta.n*log(Delta.L2/(s_temp-Delta.L2));
            
            F.function = &f2;
            s_temp = Delta.L2-eps;
            F.params = &s_temp;
            
            gsl_integration_qawc(&F,Delta.s_th,Delta.L2,s.re,eps1,eps2,M,w,&result1,&error);
            
            F.function = &high_energy_reminder;
            F.params = &s_temp;
            
            gsl_integration_qags(&F,0.,1.,eps1,eps2,M,w,&result2,&error);
            
            result1 = ((result1-result2)*s.re+(log(Delta.s_th/(s.re-Delta.s_th))+log((Delta.L2-s_temp)/Delta.L2))*delta(s.re))/M_PI+Delta.n*log(Delta.L2/(Delta.L2-s_temp));
            
            result = 0.5*(result+result1);
        }
        
        else if (s.re>Delta.L2+eps) {
            F.function = &f1;
            F.params = &s.re;
            
            gsl_integration_qags(&F,Delta.s_th,Delta.L2,eps1,eps2,M,w,&result,&error);
            
            F.function = &high_energy_reminder_cv;
            s_temp = 1.-Delta.L2/s.re;
            F.params = &s_temp;
            
            gsl_integration_qawc(&F,0.,1.,s_temp,eps1,eps2,M,w,&result1,&error);
            
            result = s.re/M_PI*(result-result1)+Delta.n*log(Delta.L2/(s.re-Delta.L2));
        }
        
        if (s.re>Delta.s_th) {
            temp.re = cos(delta(s.re))*exp(result);
            temp.im = sin(delta(s.re))*exp(result);
        }
        else {
            temp.re = exp(result);
            temp.im = 0.;
        }
        
        if (s.im<0.) {
            temp.im *= -1.;
        }
    }
    
    else {
        double f3 (double z, void * params) {
            complex x = *(complex *) params;
            double fct = delta(z)*(1.0-x.re/z)/(pow(z-x.re,2.0)+x.im*x.im);
            return fct;
        }
        
        double f4 (double z, void * params) {
            complex x = *(complex *) params;
            double fct = delta(z)*x.im/(z*(pow(z-x.re,2.0)+x.im*x.im));
            return fct;
        }
        
        double high_energy_cmp_re(double z, void * params) {
            complex x = *(complex *)params;
            double fct = Delta.a*pow(1.-z,Delta.m)*(Delta.L2-x.re*(1.-z))/((Delta.b*pow(1.-z,Delta.m)+1.)*(pow(Delta.L2-x.re*(1.-z),2.)+pow(x.im*(1.-z),2.)));
            return fct;
        }
        
        double high_energy_cmp_im(double z, void * params) {
            complex x = *(complex *)params;
            double fct = Delta.a*pow(1.-z,Delta.m)*(x.im*(1.-z))/((Delta.b*pow(1.-z,Delta.m)+1.)*(pow(Delta.L2-x.re*(1.-z),2.)+pow(x.im*(1.-z),2.)));
            return fct;
        }
        
        F.function = &f3;
        F.params = &s;
        
        gsl_integration_qags(&F,Delta.s_th,Delta.L2,eps1,eps2,M,w,&A,&error);
        
        F.function = &f4;
        F.params = &s;
        
        gsl_integration_qags(&F,Delta.s_th,Delta.L2,eps1,eps2,M,w,&B,&error);
        
        F.function = &high_energy_cmp_re;
        F.params = &s;
        
        gsl_integration_qags(&F,0.,1.,eps1,eps2,M,w,&C,&error);
        
        F.function = &high_energy_cmp_im;
        F.params = &s;
        
        gsl_integration_qags(&F,0.,1.,eps1,eps2,M,w,&D,&error);
        
        r = sqrt(pow(Delta.L2-s.re,2.0)+pow(s.im,2.0))/Delta.L2;
        phi = atan(-s.im/(Delta.L2-s.re));
        
        result = exp((s.re*(A-C)-s.im*(B-D))/M_PI-Delta.n*log(r));
        result1 = (s.im*(A-C)+s.re*(B-D))/M_PI-Delta.n*phi;
        
        temp.re = result*cos(result1);
        temp.im = result*sin(result1);
        
        if (s.re>Delta.L2) {
            temp.re *= -1.;
            temp.im *= -1.;
        }
    }
    
    gsl_integration_workspace_free(w);
    
    return temp;
}

//omnes function for given s in [s0,L2] above the cut
void omnes_cv_plus(complex *omnes, double (*delta)(double), delta_params Delta, double *s, int N) {
    int i;
    complex s_cmp;
    
    s_cmp.im = 1.e-10;
    
    for (i=0; i<N; i++) {
        s_cmp.re = s[i];
        
        omnes[i] = Omnes_function(delta,Delta,s_cmp);
        
        if (isnan(omnes[i].re)==1) {
            printf("FATAL ERROR: In function omnes_cv_plus, omnes[%d].re=nan, s=%.10e\n",i,s[i]);
            exit(1);
        }
        else if (isnan(omnes[i].im)==1) {
            printf("FATAL ERROR: In function omnes_cv_plus, omnes[%d].im=nan, s=%.10e\n",i,s[i]);
            exit(1);
        }
    }
}

//omnes function for given s in [s0,L2] below the cut
void omnes_cv_minus(complex *omnes, double (*delta)(double), delta_params Delta, double *s, int N) {
    int i;
    complex s_cmp;
    
    s_cmp.im = -1.e-10;
    
    for (i=0; i<N; i++) {
        s_cmp.re = s[i];
        
        omnes[i] = Omnes_function(delta,Delta,s_cmp);
        
        if (isnan(omnes[i].re)==1) {
            printf("FATAL ERROR: In function omnes_cv_minus, omnes[%d].re=nan, s=%.10e\n",i,s[i]);
            exit(1);
        }
        else if (isnan(omnes[i].im)==1) {
            printf("FATAL ERROR: In function omnes_cv_minus, omnes[%d].im=nan, s=%.10e\n",i,s[i]);
            exit(1);
        }
    }
}

//omnes function for given complex s in complex path for integration region III
void omnes_complex(complex *omnes, double (*delta)(double), delta_params Delta, complex *s, int N) {
    int i;
    
    for (i=0; i<N; i++) {
        omnes[i] = Omnes_function(delta,Delta,s[i]);
        
        if (isnan(omnes[i].re)==1) {
            printf("FATAL ERROR: In function omnes_complex, omnes[%d].re=nan, s.re=%.10e s.im=%.10e\n",i,s[i].re,s[i].im);
            exit(1);
        }
        else if (isnan(omnes[i].im)==1) {
            printf("FATAL ERROR: In function omnes_complex, omnes[%d].im=nan, s.re=%.10e s.im=%.10e\n",i,s[i].re,s[i].im);
            exit(1);
        }
    }
}

//omnes function for given s in [si,s0]
void omnes_below_cut(complex *omnes, double (*delta)(double), delta_params Delta, double *s, int N) {
    int i;
    complex s_cmp;
    
    s_cmp.im = 0.;
    
    for (i=0; i<N; i++) {
        s_cmp.re = s[i];
        
        omnes[i] = Omnes_function(delta,Delta,s_cmp);
        
        if (isnan(omnes[i].re)==1) {
            printf("FATAL ERROR: In function omnes_below_cut, omnes[%d].re=nan, s=%.10e\n",i,s[i]);
            exit(1);
        }
        else if (isnan(omnes[i].im)==1) {
            printf("FATAL ERROR: In function omnes_below_cut, omnes[%d].im=nan, s=%.10e\n",i,s[i]);
            exit(1);
        }
    }
}

//calculating omnes function in all four regions
void build_omnes(complex **omnes, double (*delta)(double), delta_params Delta, double *s_cv_plus, double *s_cv_minus, double *s_below_cut, complex **s_cmp, int *N) {
    omnes_cv_plus(omnes[0],delta,Delta,s_cv_plus,N[0]);
    omnes_cv_minus(omnes[1],delta,Delta,s_cv_minus,N[1]);
    omnes_below_cut(omnes[2],delta,Delta,s_below_cut,N[2]);
    omnes_complex(omnes[3],delta,Delta,s_cmp[0],N[3]);
    omnes_complex(omnes[4],delta,Delta,s_cmp[1],N[3]);
}

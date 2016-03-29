#include "Basic.h"
#include "AngularAverages.h"
#include "Omnes.h"

// delta - phaseshift delta(s) with s0 <= s <= L2
// s - Mandelstam s
// s0 - threshold
// L2 - cutoff Lambda^2
// delta_L2 - value of phaseshift delta(L2)

// TODO: class (?) for boundary and splines and just give a reference to these before
complex Omnes_function(gsl_spline *delta, complex s, double s0, double L2, double delta_L2) {
    int M = 1500;
    complex temp;
    double A,B,r,phi,result,result1,error,eps,eps1,eps2;
    
    gsl_interp_accel *acc = gsl_interp_accel_alloc();
    gsl_integration_workspace *w = gsl_integration_workspace_alloc(M);
    gsl_function F;
    
    eps = 1.0e-10;
    eps1 = 1.0e-07;
    eps2 = 1.0e-07;
    
    if (s.im<eps && s.im>-eps) {
        double f1 (double z, void * params) {
            double x = *(double *) params;
            double fct = gsl_spline_eval(delta,z,acc)/(z*(z-x));
            return fct;
        }
        
        double f2 (double z, void * params) {
            double x = *(double *) params;
            double fct = (gsl_spline_eval(delta,z,acc)-gsl_spline_eval(delta,x,acc))/z;
            return fct;
        }
        
        if (s.re<s0-eps) {
            F.function = &f1;
            F.params = &s.re;
            
            gsl_integration_qags(&F,s0,L2,eps1,eps2,M,w,&result,&error);
            
            result = (result*s.re+delta_L2*log(L2/(L2-s.re)))/M_PI;
        }
        
        else if (s.re>s0-eps && s.re<s0+eps) {
            F.function = &f1;
            s.re = s0-eps;
            F.params = &s.re;
            
            gsl_integration_qags(&F,s0,L2,eps1,eps2,M,w,&result,&error);
            
            result = (result*s.re+delta_L2*log(L2/(L2-s.re)))/M_PI;
            
            F.function = &f2;
            s.re = s0+eps;
            F.params = &s.re;
            
            gsl_integration_qawc(&F,s0,L2,s.re,eps1,eps2,M,w,&result1,&error);
            
            result1 = (result1*s.re+(log(s0/(s.re-s0))+log((L2-s.re)/L2))*gsl_spline_eval(delta,s.re,acc)+delta_L2*log(L2/(L2-s.re)))/M_PI;
            
            result = 0.5*(result+result1);
        }
        
        else if (s.re>s0+eps && s.re<L2-eps) {
            F.function = &f2;
            F.params = &s.re;
            
            gsl_integration_qawc(&F,s0,L2,s.re,eps1,eps2,M,w,&result,&error);
            
            result = (result*s.re+(log(s0/(s.re-s0))+log((L2-s.re)/L2))*gsl_spline_eval(delta,s.re,acc)+delta_L2*log(L2/(L2-s.re)))/M_PI;
        }
        
        else if (s.re>L2-eps && s.re<L2+eps) {
            F.function = &f1;
            s.re = L2+eps;
            F.params = &s.re;
            
            gsl_integration_qags(&F,s0,L2,eps1,eps2,M,w,&result,&error);
            
            result = (result*s.re+delta_L2*log(L2/(L2-s.re)))/M_PI;
            
            F.function = &f2;
            s.re = L2-eps;
            F.params = &s.re;
            
            gsl_integration_qawc(&F,s0,L2,s.re,eps1,eps2,M,w,&result1,&error);
            
            result1 = (result1*s.re+(log(s0/(s.re-s0))+log((L2-s.re)/L2))*gsl_spline_eval(delta,s.re,acc)+delta_L2*log(L2/(L2-s.re)))/M_PI;
            
            result = 0.5*(result+result1);
        }
        
        else if (s.re>L2+eps) {
            F.function = &f1;
            F.params = &s.re;
            
            gsl_integration_qags(&F,s0,L2,eps1,eps2,M,w,&result,&error);
            
            result = (result*s.re+delta_L2*log(L2/(s.re-L2)))/M_PI;
        }
        
        if (s.re>s0) {
            temp.re = cos(gsl_spline_eval(delta,s.re,acc))*exp(result);
            temp.im = sin(gsl_spline_eval(delta,s.re,acc))*exp(result);
        }
        else {
            temp.re = exp(result);
            temp.im = 0.;
        }
        
        if (s.im<0.0) {
            temp.im *= -1.0;
        }
    }
    
    else {
        double f3 (double z, void * params) {
            complex x = *(complex *) params;
            double fct = gsl_spline_eval(delta,z,acc)*(1.0-x.re/z)/(pow(z-x.re,2.0)+x.im*x.im);
            return fct;
        }
        
        double f4 (double z, void * params) {
            complex x = *(complex *) params;
            double fct = gsl_spline_eval(delta,z,acc)*x.im/(z*(pow(z-x.re,2.0)+x.im*x.im));
            return fct;
        }
        
        F.function = &f3;
        F.params = &s;
        
        gsl_integration_qags(&F,s0,L2,eps1,eps2,M,w,&A,&error);
        
        F.function = &f4;
        F.params = &s;
        
        gsl_integration_qags(&F,s0,L2,eps1,eps2,M,w,&B,&error);
        
        r = sqrt(pow(L2-s.re,2.0)+pow(s.im,2.0))/L2;
        phi = atan(-s.im/(L2-s.re));
        
        result = exp((s.re*A-s.im*B-delta_L2*log(r))/M_PI);
        result1 = (s.im*A+s.re*B-delta_L2*phi)/M_PI;
        
        temp.re = result*cos(result1);
        temp.im = result*sin(result1);
    }
    
    gsl_interp_accel_free(acc);
    gsl_integration_workspace_free(w);
    
    return temp;
}

//abs(omnes) function for given s in [s0,L2] above the cut for M_hat integration
void abs_omnes_cv_plus(double *abs_omnes, gsl_spline *delta, double s0, double L2, double delta_L2, double *s, int N) {
    int i,M;
    double error,eps,temp,pts[3];
    
    M  = 1500;
    eps = 1.0e-10;
    
    gsl_interp_accel *acc = gsl_interp_accel_alloc();
    gsl_integration_workspace *w = gsl_integration_workspace_alloc(M);
    gsl_function F;
    
    double f (double z, void * params) {
        double x = *(double *) params;
        double fct = (gsl_spline_eval(delta,z,acc)-gsl_spline_eval(delta,x,acc))/(z*(z-x));
        return fct;
    }
    
    F.function = &f;
    
    pts[0] = s0;
    pts[2] = L2;
    
    for (i=0; i<N; i++) {
        F.params = &s[i];
        
        pts[1] = s[i];
        
        gsl_integration_qagp(&F,pts,3,eps,eps,M,w,&temp,&error);
        
        abs_omnes[i] = exp((temp*s[i]+(log(s0/(s[i]-s0))+log((L2-s[i])/L2))*gsl_spline_eval(delta,s[i],acc)+delta_L2*log(L2/(L2-s[i])))/M_PI);
        
        if (isnan(abs_omnes[i])==1) {
            printf("FATAL ERROR: In function abs_omnes_cv_plus, abs_omnes[%d]=nan, s=%.10e\n",i,s[i]);
            exit(1);
        }
    }
    
    gsl_interp_accel_free(acc);
    gsl_integration_workspace_free(w);
}

//omnes function for given s in [s0,L2] above the cut
void omnes_cv_plus(complex *omnes, double *abs_omnes, gsl_spline *delta, double s0, double L2, double delta_L2, double *s, int N) {
    int i,M;
    double error,eps,temp,pts[3];
    
    M  = 1500;
    eps = 1.0e-10;
    
    gsl_interp_accel *acc = gsl_interp_accel_alloc();
    gsl_integration_workspace *w = gsl_integration_workspace_alloc(M);
    gsl_function F;
    
    double f (double z, void * params) {
        double x = *(double *) params;
        double fct = (gsl_spline_eval(delta,z,acc)-gsl_spline_eval(delta,x,acc))/(z*(z-x));
        if (isnan(fct)==1) {
            printf("%.10e %.10e\n",z,x);
        }
        return fct;
    }
    
    F.function = &f;
    
    pts[0] = s0;
    pts[2] = L2;
    
    for (i=0; i<N; i++) {
        F.params = &s[i];
        
        pts[1] = s[i];
        
        gsl_integration_qagp(&F,pts,3,eps,eps,M,w,&temp,&error);
        
        temp = exp((temp*s[i]+(log(s0/(s[i]-s0))+log((L2-s[i])/L2))*gsl_spline_eval(delta,s[i],acc)+delta_L2*log(L2/(L2-s[i])))/M_PI);
        omnes[i].re = cos(gsl_spline_eval(delta,s[i],acc))*temp;
        omnes[i].im = sin(gsl_spline_eval(delta,s[i],acc))*temp;
        
        abs_omnes[i] = temp;
        
        if (isnan(abs_omnes[i])==1) {
            printf("FATAL ERROR: In function omnes_cv_plus, abs_omnes[%d]=nan, s=%.10e\n",i,s[i]);
            exit(1);
        }
        else if (isnan(omnes[i].re)==1) {
            printf("FATAL ERROR: In function omnes_cv_plus, omnes[%d].re=nan, s=%.10e\n",i,s[i]);
            exit(1);
        }
        else if (isnan(omnes[i].im)==1) {
            printf("FATAL ERROR: In function omnes_cv_plus, omnes[%d].im=nan, s=%.10e\n",i,s[i]);
            exit(1);
        }
    }
    
    gsl_interp_accel_free(acc);
    gsl_integration_workspace_free(w);
}

//omnes function for given s in [s0,L2] below the cut
void omnes_cv_minus(complex *omnes, gsl_spline *delta, double s0, double L2, double delta_L2, double *s, int N) {
    int i,M;
    double error,eps,temp,pts[3];
    
    M  = 1500;
    eps = 1.0e-10;
    
    gsl_interp_accel *acc = gsl_interp_accel_alloc();
    gsl_integration_workspace *w = gsl_integration_workspace_alloc(M);
    gsl_function F;
    
    double f (double z, void * params) {
        double x = *(double *) params;
        double fct = (gsl_spline_eval(delta,z,acc)-gsl_spline_eval(delta,x,acc))/(z*(z-x));
        return fct;
    }
    
    F.function = &f;
    
    pts[0] = s0;
    pts[2] = L2;
    
    for (i=0; i<N; i++) {
        F.params = &s[i];
        
        pts[1] = s[i];
        
        gsl_integration_qagp(&F,pts,3,eps,eps,M,w,&temp,&error);
        
        temp = exp((temp*s[i]+(log(s0/(s[i]-s0))+log((L2-s[i])/L2))*gsl_spline_eval(delta,s[i],acc)+delta_L2*log(L2/(L2-s[i])))/M_PI);
        omnes[i].re = cos(gsl_spline_eval(delta,s[i],acc))*temp;
        omnes[i].im = -sin(gsl_spline_eval(delta,s[i],acc))*temp;
        
        if (isnan(omnes[i].re)==1) {
            printf("FATAL ERROR: In function omnes_cv_minus, omnes[%d].re=nan, s=%.10e\n",i,s[i]);
            exit(1);
        }
        else if (isnan(omnes[i].im)==1) {
            printf("FATAL ERROR: In function omnes_cv_minus, omnes[%d].im=nan, s=%.10e\n",i,s[i]);
            exit(1);
        }
    }
    
    gsl_interp_accel_free(acc);
    gsl_integration_workspace_free(w);
}

//omnes function for given complex s in complex path for integration region III
void omnes_complex(complex *omnes, gsl_spline *delta, double s0, double L2, double delta_L2, complex *s, int N) {
    int i,M;
    double error,eps,A,B,r,gamma;
    
    M  = 1500;
    eps = 1.0e-10;
    
    gsl_interp_accel *acc = gsl_interp_accel_alloc();
    gsl_integration_workspace *w = gsl_integration_workspace_alloc(M);
    gsl_function F;
    
    double f_re (double z, void * params) {
        complex x = *(complex *) params;
        double fct = gsl_spline_eval(delta,z,acc)*(1.0-x.re/z)/(pow(z-x.re,2.0)+x.im*x.im);
        return fct;
    }
    
    double f_im (double z, void * params) {
        complex x = *(complex *) params;
        double fct = gsl_spline_eval(delta,z,acc)*x.im/(z*(pow(z-x.re,2.0)+x.im*x.im));
        return fct;
    }
    
    for (i=0; i<N; i++) {
        F.params = &s[i];
        
        F.function = &f_re;
        gsl_integration_qags(&F,s0,L2,eps,eps,M,w,&A,&error);
        
        F.function = &f_im;
        gsl_integration_qags(&F,s0,L2,eps,eps,M,w,&B,&error);
        
        r = sqrt(pow(L2-s[i].re,2.0)+pow(s[i].im,2.0))/L2;
        gamma = atan(-s[i].im/(L2-s[i].re));
        
        r = exp((s[i].re*A-s[i].im*B-delta_L2*log(r))/M_PI);
        gamma = (s[i].im*A+s[i].re*B-delta_L2*gamma)/M_PI;
        
        omnes[i].re = r*cos(gamma);
        omnes[i].im = r*sin(gamma);
        
        if (isnan(omnes[i].re)==1) {
            printf("FATAL ERROR: In function omnes_complex, omnes[%d].re=nan, s.re=%.10e s.im=%.10e\n",i,s[i].re,s[i].im);
            exit(1);
        }
        else if (isnan(omnes[i].im)==1) {
            printf("FATAL ERROR: In function omnes_complex, omnes[%d].im=nan, s.re=%.10e s.im=%.10e\n",i,s[i].re,s[i].im);
            exit(1);
        }
    }
    
    gsl_interp_accel_free(acc);
    gsl_integration_workspace_free(w);
}

//omnes function for given s in [si,s0]
void omnes_below_cut(complex *omnes, gsl_spline *delta, double s0, double L2, double delta_L2, double *s, int N) {
    int i,M;
    double error,eps,temp;
    
    M  = 1500;
    eps = 1.0e-10;
    
    gsl_interp_accel *acc = gsl_interp_accel_alloc();
    gsl_integration_workspace *w = gsl_integration_workspace_alloc(M);
    gsl_function F;
    
    double f (double z, void * params) {
        double x = *(double *) params;
        double fct = gsl_spline_eval(delta,z,acc)/(z*(z-x));
        return fct;
    }
    
    F.function = &f;
    
    for (i=0; i<N; i++) {
        F.params = &s[i];
        
        gsl_integration_qags(&F,s0,L2,eps,eps,M,w,&temp,&error);
        
        temp = exp((temp*s[i]+delta_L2*log(L2/(L2-s[i])))/M_PI);
        omnes[i].re = temp;
        omnes[i].im = 0.0;
        
        if (isnan(omnes[i].re)==1) {
            printf("FATAL ERROR: In function omnes_below_cut, omnes[%d].re=nan, s=%.10e\n",i,s[i]);
            exit(1);
        }
        else if (isnan(omnes[i].im)==1) {
            printf("FATAL ERROR: In function omnes_below_cut, omnes[%d].im=nan, s=%.10e\n",i,s[i]);
            exit(1);
        }
    }
    
    gsl_interp_accel_free(acc);
    gsl_integration_workspace_free(w);
}

//calculating omnes function in all four regions
void build_omnes(complex **omnes, double *abs_omnes, gsl_spline *delta, double s0, double L2, double delta_L2, double *s_cv_plus, double *s_cv_minus, double *s_below_cut, complex *s_complex, int *N) {
    omnes_cv_plus(omnes[0],abs_omnes,delta,s0,L2,delta_L2,s_cv_plus,N[0]);
    omnes_cv_minus(omnes[1],delta,s0,L2,delta_L2,s_cv_minus,N[1]);
    omnes_complex(omnes[2],delta,s0,L2,delta_L2,s_complex,N[2]);
    omnes_below_cut(omnes[3],delta,s0,L2,delta_L2,s_below_cut,N[3]);
}

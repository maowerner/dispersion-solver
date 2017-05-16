#include "stdio.h"
#include <math.h>
#include <gsl/gsl_integration.h>

#define a 34.1389795918
#define b 61.5104081633
#define c 46.8246938776
#define const1 7.84285714286
#define const2 -5.84285714286

typedef  struct {
    double re, im;
}complex;

double integration_III_radius_phi(double phi) {
    double result,z,eta,eps;
    
    eps = 1.0e-05;
    
    z = cos(phi);
    
    if (z==-1.0) {
        result = sqrt(a);
    }
    else if (z>-1.0 && z<-eps) {
        eta = 1.0/3.0*(acos(54.0*pow(c-1.0,2.0)/pow(c+3.0,3.0)*pow(z,2.0)-1.0));
        result = (c+3.0)/(6.0*z)*(1.0-2.0*cos(eta));
    }
    else if (z>=-eps && z<=eps) {
        result = (c-1.0)/sqrt(c+3.0);
    }
    else if (z>eps && z<1.0) {
        eta = 1.0/3.0*(acos(54.0*pow(c-1.0,2.0)/pow(c+3.0,3.0)*pow(z,2.0)-1.0)+M_PI);
        result = (c+3.0)/(6.0*z)*(1.0+2.0*cos(eta));
    }
    else if (z==1.0) {
        result = sqrt(b);
    }
    
    return result;
}

complex integration_III_gamma_phi(double phi) {
    double r;
    complex result;
    
    r = integration_III_radius_phi(phi);
    
    result.re = r*cos(phi);
    result.im = r*sin(phi);
    
    return result;
}

double Q12_re_int(double z, void * params) {
    double s = *(double *) params;
    double f;
    
    if (s<a) {
        if (z<a) {
            f = 1.0/(sqrt(a-z)*(z-s));
        }
        else {
            f = 0.0;
        }
    }
    else {
        if (z<a) {
            f = 0.0;
        }
        else {
            f = 1.0/(sqrt(z-a)*(z-s));
        }
    }
    
    return f;
}

double Q12_im_int(double z, void * params) {
    double s = *(double *) params;
    double f;
    
    if (s<a) {
        if (z<a) {
            f = 0.0;
        }
        else {
            f = 1.0/(sqrt(z-a)*(z-s)); //-
        }
    }
    else {
        if (z<a) {
            f = 1.0/(sqrt(a-z)*(z-s));
        }
        else {
            f = 0.0;
        }
    }
    
    return f;
}

complex Q12_real(double s, double x, double y) {
    int w_size = 1000;
    double error,abserr,relerr,pts[3];
    complex result;
    
    abserr = 1.0e-10;
    relerr = 1.0e-10;
    
    pts[0] = x;
    pts[1] = a;
    pts[2] = y;
    
    gsl_integration_workspace *w = gsl_integration_workspace_alloc(w_size);
    
    gsl_function F;
    F.params = &s;
    
    F.function = &Q12_re_int;
    gsl_integration_qagp(&F,pts,3,abserr,relerr,w_size,w,&result.re,&error);
    
    F.function = &Q12_im_int;
    gsl_integration_qagp(&F,pts,3,abserr,relerr,w_size,w,&result.im,&error);
    
    gsl_integration_workspace_free(w);
    
    return result;
}

double R12_int(double z, void * params) {
    double s = *(double *) params;
    double f;
    
    if (s<a) {
        f = 1.0/sqrt(a-z);
    }
    else {
        f = 1.0/sqrt(z-a);
    }
    
    return f;
}

double R12_real(double s, double x, double y) {
    int w_size = 1000;
    double result,error,abserr,relerr;
    
    abserr = 1.0e-10;
    relerr = 1.0e-10;
    
    gsl_integration_workspace *w = gsl_integration_workspace_alloc(w_size);
    
    gsl_function F;
    F.params = &s;
    F.function = &R12_int;
    gsl_integration_qawc(&F,x,y,s,abserr,relerr,w_size,w,&result,&error);
    
    gsl_integration_workspace_free(w);
    
    return result;
}

double R32_int(double z, void * params) {
    double s = *(double *) params;
    double f;
    
    if (s<a) {
        f = 1.0/pow(a-z,1.5);
    }
    else {
        f = 1.0/pow(z-a,1.5);
    }
    
    return f;
}

double R32_real(double s, double x, double y) {
    int w_size = 1000;
    double result,error,abserr,relerr;
    
    abserr = 1.0e-10;
    relerr = 1.0e-10;
    
    gsl_integration_workspace *w = gsl_integration_workspace_alloc(w_size);
    
    gsl_function F;
    F.params = &s;
    F.function = &R32_int;
    gsl_integration_qawc(&F,x,y,s,abserr,relerr,w_size,w,&result,&error);
    
    gsl_integration_workspace_free(w);
    
    return result;
}

double Q12_re_int_cmp(double z, void * params) {
    complex s = *(complex *) params;
    double x,y,f;
    
    x = z-s.re;
    y = x*x+s.im*s.im;
    
    if (z<a) {
        f = 0.0;//x/(pow(a-z,1.5)*y);
    }
    else {
        f = x/(pow(z-a,1.5)*y);
    }
    
    return f;
}

double Q12_im_int_cmp(double z, void * params) {
    complex s = *(complex *) params;
    double x,y,f;
    
    x = z-s.re;
    y = x*x+s.im*s.im;
    
    if (z<a) {
        f = 0.0;//s.im/(pow(a-z,1.5)*y);
    }
    else {
        f = s.im/(pow(z-a,1.5)*y);
    }
    
    return f;
}

complex Q12_cmp(complex s, double x, double y) {
    int w_size = 1000;
    double error,abserr,relerr,pts[3];
    complex result;
    
    abserr = 1.0e-10;
    relerr = 1.0e-10;
    
    pts[0] = x;
    pts[1] = a;
    pts[2] = y;
    
    gsl_integration_workspace *w = gsl_integration_workspace_alloc(w_size);
    
    gsl_function F;
    F.params = &s;
    
    F.function = &Q12_re_int_cmp;
    gsl_integration_qagp(&F,pts,3,abserr,relerr,w_size,w,&result.re,&error);
    
    F.function = &Q12_im_int_cmp;
    gsl_integration_qagp(&F,pts,3,abserr,relerr,w_size,w,&result.im,&error);
    
    gsl_integration_workspace_free(w);
    
    return result;
}

double Q32_re_int(double z, void * params) {
    double s = *(double *) params;
    double f;
    
    if (s<a) {
        if (z<a) {
            f = 1.0/(pow(a-z,1.5)*(z-s));
        }
        else {
            f = 0.0;
        }
    }
    else {
        if (z<a) {
            f = 0.0;
        }
        else {
            f = 1.0/(pow(z-a,1.5)*(z-s));
        }
    }
    
    return f;
}

double Q32_im_int(double z, void * params) {
    double s = *(double *) params;
    double f;
    
    if (s<a) {
        if (z<a) {
            f = 0.0;
        }
        else {
            f = 1.0/(pow(z-a,1.5)*(z-s));
        }
    }
    else {
        if (z<a) {
            f = 1.0/(pow(a-z,1.5)*(z-s));//-
        }
        else {
            f = 0.0;
        }
    }
    
    return f;
}

complex Q32_real(double s, double x, double y) {
    int w_size = 1000;
    double error,abserr,relerr,pts[3];
    complex result;
    
    abserr = 1.0e-5;
    relerr = 1.0e-5;
    
    pts[0] = x;
    pts[1] = a;
    pts[2] = y;
    
    gsl_integration_workspace *w = gsl_integration_workspace_alloc(w_size);
    
    gsl_function F;
    F.params = &s;
    
    F.function = &Q32_re_int;
    gsl_integration_qagp(&F,pts,3,abserr,relerr,w_size,w,&result.re,&error);
    
    F.function = &Q32_im_int;
    gsl_integration_qagp(&F,pts,3,abserr,relerr,w_size,w,&result.im,&error);
    
    gsl_integration_workspace_free(w);
    
    return result;
}

int main() {
    int i,N;
    double s,ds,si,sf,s0,L2,eps,R12,R32;
    complex Q12,Q32,gamma;
    
    gsl_set_error_handler_off();
    
    s0 = 4.0;
    L2 = 150.0;
    
    eps = 1.0e-1;
    
    N = 10000;
    si = 0.0;
    sf = 2.0*M_PI;
    ds = (sf-si)/(double)(N-1);
    for (i=0; i<N; i++) {
        s = si+ds*(double)i;
        gamma = integration_III_gamma_phi(s);
        Q12 = Q12_cmp(gamma,s0,L2);
        //Q32 = Q32_real(gamma,s0,L2);
        printf("%.6e %.6e %.6e %.6e %.6e %.6e %.6e\n",s,gamma.re,gamma.im,Q12.re,Q12.im);
    }
    
    /*gamma.re = const1;
    gamma.im = 1.0e-5;
    Q12 = Q12_cmp(gamma,s0,L2);
    Q32 = Q12_real(const1,0.5*(a+const1),L2);
    Q32.re += R12_real(const1,s0,0.5*(a+const1));
    Q32.im += M_PI/sqrt(a-const1);
    printf("%.6e %.6e %.6e %.6e %.6e %.6e\n",Q12.re,Q12.im,Q32.re,Q32.im,fabs(Q12.re-Q32.re),fabs(Q12.im-Q32.im));
    
    gamma.re = const2;
    gamma.im = 0.0;
    Q12 = Q12_cmp(gamma,s0,L2);
    Q32 = Q12_real(const2,s0,L2);
    printf("%.6e %.6e %.6e %.6e %.6e %.6e\n",Q12.re,Q12.im,Q32.re,Q32.im,fabs(Q12.re-Q32.re),fabs(Q12.im-Q32.im));*/
    
//    N = 10000;
//    si = -100.0;
//    sf = 4.0-eps;
//    ds = (sf-si)/(double)(N-1);
//    for (i=0; i<N; i++) {
//        s = si+ds*(double)i;
//        Q12 = Q12_real(s,s0,L2);
//        Q32 = Q32_real(s,s0,L2);
//        printf("%.6e %.6e %.6e %.6e %.6e\n",s,Q12.re,Q12.im,Q32.re,Q32.im);
//    }
//    
//    N = 10000;
//    si = 4.0+eps;
//    sf = a-eps;
//    ds = (sf-si)/(double)(N-1);
//    for (i=0; i<N; i++) {
//        s = si+ds*(double)i;
//        Q12 = Q12_real(s,0.5*(a+s),L2);
//        R12 = R12_real(s,s0,0.5*(a+s));
//        Q32 = Q32_real(s,0.5*(a+s),L2);
//        R32 = R32_real(s,s0,0.5*(a+s));
//        printf("%.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e\n",s,Q12.re,Q12.im,Q32.re,Q32.im,R12,M_PI/sqrt(a-s),R32,M_PI/pow(a-s,1.5));
//    }
//    
//    N = 10000;
//    si = a+eps;
//    sf = L2-eps;
//    ds = (sf-si)/(double)(N-1);
//    for (i=0; i<N; i++) {
//        s = si+ds*(double)i;
//        Q12 = Q12_real(s,s0,0.5*(a+s));
//        R12 = R12_real(s,0.5*(a+s),L2);
//        Q32 = Q32_real(s,s0,0.5*(a+s));
//        R32 = R32_real(s,0.5*(a+s),L2);
//        printf("%.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e\n",s,Q12.re,Q12.im,Q32.re,Q32.im,R12,M_PI/sqrt(s-a),R32,M_PI/pow(s-a,1.5));
//    }
    
    return 0;
}

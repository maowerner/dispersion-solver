#include "Basic.h"
#include "InhomogenitySqrtCubed.h"

//
void inhomogenity_sqrt_cubed_cv_plus(complex *M_inhom, complex_spline *M_hat, complex *f, complex *g, double *s, complex *R32, complex *Q12, complex *Q32, double s0, double L2, double cut, int N, int n) {
    int i,N_int,Nm,Np;
    double a,am,ap,p,eps,error,interpol_re[4],interpol_im[4],s3,s2,pts_cv[2][3],pts_a[2][3],log_const;
    complex result,M,dM,M_const;
    
    a = METAP_M_MPION_SQUARED;
    
    N_int = 1500;
    eps = 1.0e-14;
    
    am = a-DELTA_A;
    ap = a+DELTA_A;
    
    gsl_interp_accel *acc0 = gsl_interp_accel_alloc();
    gsl_interp_accel *acc1 = gsl_interp_accel_alloc();
    gsl_interp_accel *acc2 = gsl_interp_accel_alloc();
    
    gsl_integration_workspace *w = gsl_integration_workspace_alloc(N_int);
    
    gsl_function F_re_cv;
    gsl_function F_im_cv;
    gsl_function F_re_sing;
    gsl_function F_im_sing;
    
    M_const.re = gsl_spline_eval(M_hat[1].re,L2,acc2)*L2/pow(L2-a,1.5);
    M_const.im = gsl_spline_eval(M_hat[1].im,L2,acc2)*L2/pow(L2-a,1.5);
    
    double f_re_cv(double z, void * params) {
        double x = *(double *) params;
        double fct;
        if (z<a) {
            fct = (gsl_spline_eval(M_hat[0].re,z,acc1)-gsl_spline_eval(M_hat[0].re,x,acc1))/(pow(a-z,1.5)*(z-x));
        }
        else {
            fct = (gsl_spline_eval(M_hat[1].re,z,acc2)-gsl_spline_eval(M_hat[1].re,x,acc2))/(pow(z-a,1.5)*(z-x));
        }
        return fct;
    }
    
    double f_im_cv(double z, void * params) {
        double x = *(double *) params;
        double fct;
        if (z<a) {
            fct = (gsl_spline_eval(M_hat[0].im,z,acc1)-gsl_spline_eval(M_hat[0].im,x,acc1))/(pow(a-z,1.5)*(z-x));
        }
        else {
            fct = (gsl_spline_eval(M_hat[1].im,z,acc2)-gsl_spline_eval(M_hat[1].im,x,acc2))/(pow(z-a,1.5)*(z-x));
        }
        return fct;
    }
    
    double f_re_sing(double z, void * params) {
        double x = *(double *) params;
        double fct;
        if (z<am) {
            fct = (gsl_spline_eval(M_hat[0].re,z,acc0)-f[0].re-f[1].re*(a-z))/(pow(a-z,1.5)*(z-x));
        }
        else if (z>=am && z<a) {
            fct = (f[2].re+f[3].re*sqrt(a-z)+f[4].re*(a-z))/(z-x);
        }
        else if (z>a && z<=ap) {
            fct = (g[2].re+g[3].re*sqrt(z-a)+g[4].re*(z-a))/(z-x);
        }
        else {
            fct = (gsl_spline_eval(M_hat[1].re,z,acc2)-g[0].re-g[1].re*(z-a))/(pow(z-a,1.5)*(z-x));
        }
        return fct;
    }
    
    double f_im_sing(double z, void * params) {
        double x = *(double *) params;
        double fct;
        if (z<am) {
            fct = (gsl_spline_eval(M_hat[0].im,z,acc1)-f[0].im-f[1].im*(a-z))/(pow(a-z,1.5)*(z-x));
        }
        else if (z>=am && z<a) {
            fct = (f[2].im+f[3].im*sqrt(a-z)+f[4].im*(a-z))/(z-x);
        }
        else if (z>a && z<=ap) {
            fct = (g[2].im+g[3].im*sqrt(z-a)+g[4].im*(z-a))/(z-x);
        }
        else {
            fct = (gsl_spline_eval(M_hat[1].im,z,acc2)-g[0].im-g[1].im*(z-a))/(pow(z-a,1.5)*(z-x));
        }
        return fct;
    }
    
    F_re_cv.function = &f_re_cv;
    F_im_cv.function = &f_im_cv;
    F_re_sing.function = &f_re_sing;
    F_im_sing.function = &f_im_sing;
    
    pts_cv[0][0] = s0;
    pts_cv[1][2] = L2;
    
    pts_a[0][0] = s0;
    pts_a[0][1] = a;
    pts_a[1][1] = a;
    pts_a[1][2] = L2;
    
    for (i=0; i<N; i++) {
        p = 0.5*(a+s[i]);
        
        F_re_cv.params = &s[i];
        F_im_cv.params = &s[i];
        
        F_re_sing.params = &s[i];
        F_im_sing.params = &s[i];
        
        if (s[i]<a-cut) {
            Nm = i;
            
            pts_cv[0][1] = s[i];
            pts_cv[0][2] = p;
            
            pts_a[1][0] = p;
            
            gsl_integration_qagp(&F_re_cv,pts_cv[0],3,eps,eps,N_int,w,&result.re,&error);
            gsl_integration_qagp(&F_im_cv,pts_cv[0],3,eps,eps,N_int,w,&result.im,&error);
            
            M_inhom[i].re = result.re;
            M_inhom[i].im = result.im;
            
            gsl_integration_qagp(&F_re_sing,pts_a[1],3,eps,eps,N_int,w,&result.re,&error);
            gsl_integration_qagp(&F_im_sing,pts_a[1],3,eps,eps,N_int,w,&result.im,&error);
            
            M_inhom[i].re += result.re;
            M_inhom[i].im += result.im;
            
            M.re = gsl_spline_eval(M_hat[0].re,s[i],acc1);
            M.im = gsl_spline_eval(M_hat[0].im,s[i],acc1);
            
            log_const = log(L2/(L2-s[i]))/s[i];
            
            M_inhom[i].re += M.re*R32[i].re-M.im*R32[i].im+M_const.re*log_const;
            M_inhom[i].im += M.im*R32[i].re+M.re*R32[i].im+M_const.im*log_const;
            
            M_inhom[i].re += f[0].re*Q32[i].re+g[0].re*Q32[i].im;
            M_inhom[i].re += f[1].re*Q12[i].re+g[1].re*Q12[i].im;
            
            M_inhom[i].im += f[0].im*Q32[i].re+g[0].im*Q32[i].im;
            M_inhom[i].im += f[1].im*Q12[i].re+g[1].im*Q12[i].im;
        }
        
        else if (s[i]>a+cut) {
            pts_cv[1][0] = p;
            pts_cv[1][1] = s[i];
            
            pts_a[0][2] = p;
            
            gsl_integration_qagp(&F_re_cv,pts_cv[1],3,eps,eps,N_int,w,&result.re,&error);
            gsl_integration_qagp(&F_im_cv,pts_cv[1],3,eps,eps,N_int,w,&result.im,&error);
            
            M_inhom[i].re = result.re;
            M_inhom[i].im = result.im;
            
            gsl_integration_qagp(&F_re_sing,pts_a[0],3,eps,eps,N_int,w,&result.re,&error);
            gsl_integration_qagp(&F_im_sing,pts_a[0],3,eps,eps,N_int,w,&result.im,&error);
            
            M_inhom[i].re += result.re;
            M_inhom[i].im += result.im;
            
            M.re = gsl_spline_eval(M_hat[1].re,s[i],acc2);
            M.im = gsl_spline_eval(M_hat[1].im,s[i],acc2);
            
            log_const = log(L2/(L2-s[i]))/s[i];
            
            M_inhom[i].re += M.re*R32[i].re-M.im*R32[i].im+M_const.re*log_const;
            M_inhom[i].im += M.im*R32[i].re+M.re*R32[i].im+M_const.im*log_const;
            
            M_inhom[i].re += g[0].re*Q32[i].re+f[0].re*Q32[i].im;
            M_inhom[i].re += g[1].re*Q12[i].re+f[1].re*Q12[i].im;
            
            M_inhom[i].im += g[0].im*Q32[i].re+f[0].im*Q32[i].im;
            M_inhom[i].im += g[1].im*Q12[i].re+f[1].im*Q12[i].im;
        }
        
        else {
            Np = i+1;
            M_inhom[i].re = 0.0/0.0;
            M_inhom[i].im = 0.0/0.0;
        }
        
        //printf("%.10e %.10e %.10e\n",s[i],M_inhom[i].re,M_inhom[i].im);
        
        M_inhom[i].re *= pow(s[i],(double)n)/M_PI;
        M_inhom[i].im *= pow(s[i],(double)n)/M_PI;
    }
    
    dM.re = 0.;
    dM.im = 0.;
    for (i=0; i<20; i++) {
        dM.re += (M_inhom[Nm-i].re-M_inhom[Nm-1-i].re)/(s[Nm-i]-s[Nm-1-i]);
        dM.im += (M_inhom[Nm-i].im-M_inhom[Nm-1-i].im)/(s[Nm-i]-s[Nm-1-i]);
    }
    dM.re /= 20.;
    dM.im /= 20.;
    
    cubic_spline_f01_df01(s[Nm],s[Np],M_inhom[Nm].re,M_inhom[Np].re,dM.re,(M_inhom[Np+1].re-M_inhom[Np].re)/(s[Np+1]-s[Np]),interpol_re);
    cubic_spline_f01_df01(s[Nm],s[Np],M_inhom[Nm].im,M_inhom[Np].im,dM.im,(M_inhom[Np+1].im-M_inhom[Np].im)/(s[Np+1]-s[Np]),interpol_im);
    
    for (i=Nm+1; i<Np; i++) {
        s2 = s[i]*s[i];
        s3 = s[i]*s2;
        M_inhom[i].re = interpol_re[0]*s3+interpol_re[1]*s2+interpol_re[2]*s[i]+interpol_re[3];
        M_inhom[i].im = interpol_im[0]*s3+interpol_im[1]*s2+interpol_im[2]*s[i]+interpol_im[3];
    }
    
    //for (i=0; i<N; i++) {
    //    printf("%.10e %.10e %.10e\n",s[i],M_inhom[i].re,M_inhom[i].im);
    //}
    
    gsl_integration_workspace_free(w);
    gsl_interp_accel_free(acc0);
    gsl_interp_accel_free(acc1);
    gsl_interp_accel_free(acc2);
}

void inhomogenity_sqrt_cubed_cv_minus(complex *M_inhom, complex_spline *M_hat, complex *f, complex *g, double *s, complex *R32, complex *Q12, complex *Q32, double s0, double L2, int N, int n) {
    int i,N_int;
    double a,am,ap,p,eps,error,pts_cv[3],pts_a[3],log_const;
    complex result,M,M_const;
    
    a = METAP_M_MPION_SQUARED;
    
    N_int = 1500;
    eps = 1.0e-12;
    
    am = a-DELTA_A;
    ap = a+DELTA_A;
    
    gsl_interp_accel *acc0 = gsl_interp_accel_alloc();
    gsl_interp_accel *acc1 = gsl_interp_accel_alloc();
    gsl_interp_accel *acc2 = gsl_interp_accel_alloc();
    
    gsl_integration_workspace *w = gsl_integration_workspace_alloc(N_int);
    
    gsl_function F_re_cv;
    gsl_function F_im_cv;
    gsl_function F_re_sing;
    gsl_function F_im_sing;
    
    M_const.re = gsl_spline_eval(M_hat[1].re,L2,acc2)*L2/pow(L2-a,1.5);
    M_const.im = gsl_spline_eval(M_hat[1].im,L2,acc2)*L2/pow(L2-a,1.5);
    
    double f_re_cv(double z, void * params) {
        double x = *(double *) params;
        double fct;
        if (z<a) {
            fct = (gsl_spline_eval(M_hat[0].re,z,acc1)-gsl_spline_eval(M_hat[0].re,x,acc1))/(pow(a-z,1.5)*(z-x));
        }
        else {
            fct = (gsl_spline_eval(M_hat[1].re,z,acc2)-gsl_spline_eval(M_hat[1].re,x,acc2))/(pow(z-a,1.5)*(z-x));
        }
        return fct;
    }
    
    double f_im_cv(double z, void * params) {
        double x = *(double *) params;
        double fct;
        if (z<a) {
            fct = (gsl_spline_eval(M_hat[0].im,z,acc1)-gsl_spline_eval(M_hat[0].im,x,acc1))/(pow(a-z,1.5)*(z-x));
        }
        else {
            fct = (gsl_spline_eval(M_hat[1].im,z,acc2)-gsl_spline_eval(M_hat[1].im,x,acc2))/(pow(z-a,1.5)*(z-x));
        }
        return fct;
    }
    
    double f_re_sing(double z, void * params) {
        double x = *(double *) params;
        double fct;
        if (z<am) {
            fct = (gsl_spline_eval(M_hat[0].re,z,acc0)-f[0].re-f[1].re*(a-z))/(pow(a-z,1.5)*(z-x));
        }
        else if (z>=am && z<a) {
            fct = (f[2].re+f[3].re*sqrt(a-z)+f[4].re*(a-z))/(z-x);
        }
        else if (z>a && z<=ap) {
            fct = (g[2].re+g[3].re*sqrt(z-a)+g[4].re*(z-a))/(z-x);
        }
        else {
            fct = (gsl_spline_eval(M_hat[1].re,z,acc2)-g[0].re-g[1].re*(z-a))/(pow(z-a,1.5)*(z-x));
        }
        return fct;
    }
    
    double f_im_sing(double z, void * params) {
        double x = *(double *) params;
        double fct;
        if (z<am) {
            fct = (gsl_spline_eval(M_hat[0].im,z,acc1)-f[0].im-f[1].im*(a-z))/(pow(a-z,1.5)*(z-x));
        }
        else if (z>=am && z<a) {
            fct = (f[2].im+f[3].im*sqrt(a-z)+f[4].im*(a-z))/(z-x);
        }
        else if (z>a && z<=ap) {
            fct = (g[2].im+g[3].im*sqrt(z-a)+g[4].im*(z-a))/(z-x);
        }
        else {
            fct = (gsl_spline_eval(M_hat[1].im,z,acc2)-g[0].im-g[1].im*(z-a))/(pow(z-a,1.5)*(z-x));
        }
        return fct;
    }
    
    F_re_cv.function = &f_re_cv;
    F_im_cv.function = &f_im_cv;
    F_re_sing.function = &f_re_sing;
    F_im_sing.function = &f_im_sing;
    
    pts_cv[0] = s0;
    
    pts_a[1] = a;
    pts_a[2] = L2;
    
    for (i=0; i<N; i++) {
        p = 0.5*(a+s[i]);
        
        pts_cv[1] = s[i];
        pts_cv[2] = p;
        
        pts_a[0] = p;
        
        F_re_cv.params = &s[i];
        F_im_cv.params = &s[i];
        
        F_re_sing.params = &s[i];
        F_im_sing.params = &s[i];
        
        gsl_integration_qagp(&F_re_cv,pts_cv,3,eps,eps,N_int,w,&result.re,&error);
        gsl_integration_qagp(&F_im_cv,pts_cv,3,eps,eps,N_int,w,&result.im,&error);
        
        M_inhom[i].re = result.re;
        M_inhom[i].im = result.im;
        
        gsl_integration_qagp(&F_re_sing,pts_a,3,eps,eps,N_int,w,&result.re,&error);
        gsl_integration_qagp(&F_im_sing,pts_a,3,eps,eps,N_int,w,&result.im,&error);
        
        log_const = log(L2/(L2-s[i]))/s[i];
        
        M_inhom[i].re += result.re+M_const.re*log_const;
        M_inhom[i].im += result.im+M_const.im*log_const;
        
        M.re = gsl_spline_eval(M_hat[0].re,s[i],acc1);
        M.im = gsl_spline_eval(M_hat[0].im,s[i],acc1);
        
        M_inhom[i].re += M.re*R32[i].re+M.im*R32[i].im;
        M_inhom[i].im += M.im*R32[i].re-M.re*R32[i].im;
        
        M_inhom[i].re += f[0].re*Q32[i].re+g[0].re*Q32[i].im;
        M_inhom[i].re += f[1].re*Q12[i].re+g[1].re*Q12[i].im;
        
        M_inhom[i].im += f[0].im*Q32[i].re+g[0].im*Q32[i].im;
        M_inhom[i].im += f[1].im*Q12[i].re+g[1].im*Q12[i].im;
        
        M_inhom[i].re *= pow(s[i],(double)n)/M_PI;
        M_inhom[i].im *= pow(s[i],(double)n)/M_PI;
    }
    
    gsl_integration_workspace_free(w);
    gsl_interp_accel_free(acc1);
    gsl_interp_accel_free(acc2);
}

void inhomogenity_sqrt_cubed_below_cut(complex *M_inhom, complex_spline *M_hat, complex *f, complex *g, double *s, complex *Q12, complex *Q32, double s0, double L2, int N, int n) {
    int i,N_int;
    double a,am,ap,eps,error,pts[3],log_const;
    complex result,M_const;
    
    a = METAP_M_MPION_SQUARED;
    
    double *M_re = (double *)malloc(N*sizeof(double));
    double *M_im = (double *)malloc(N*sizeof(double));
    
    N_int = 1500;
    eps = 1.0e-12;
    
    am = a-DELTA_A;
    ap = a+DELTA_A;
    
    gsl_interp_accel *acc1 = gsl_interp_accel_alloc();
    gsl_interp_accel *acc2 = gsl_interp_accel_alloc();
    
    gsl_integration_workspace *w = gsl_integration_workspace_alloc(N_int);
    
    gsl_function F_re_sing;
    gsl_function F_im_sing;
    
    M_const.re = gsl_spline_eval(M_hat[1].re,L2,acc2)*L2/pow(L2-a,1.5);
    M_const.im = gsl_spline_eval(M_hat[1].im,L2,acc2)*L2/pow(L2-a,1.5);
    
    double f_re_sing(double z, void * params) {
        double x = *(double *) params;
        double fct;
        if (z<am) {
            fct = (gsl_spline_eval(M_hat[0].re,z,acc1)-f[0].re-f[1].re*(a-z))/(pow(a-z,1.5)*(z-x));
        }
        else if (z>=am && z<a) {
            fct = (f[2].re+f[3].re*sqrt(a-z)+f[4].re*(a-z))/(z-x);
        }
        else if (z>a && z<=ap) {
            fct = (g[2].re+g[3].re*sqrt(z-a)+g[4].re*(z-a))/(z-x);
        }
        else {
            fct = (gsl_spline_eval(M_hat[1].re,z,acc2)-g[0].re-g[1].re*(z-a))/(pow(z-a,1.5)*(z-x));
        }
        return fct;
    }
    
    double f_im_sing(double z, void * params) {
        double x = *(double *) params;
        double fct;
        if (z<am) {
            fct = (gsl_spline_eval(M_hat[0].im,z,acc1)-f[0].im-f[1].im*(a-z))/(pow(a-z,1.5)*(z-x));
        }
        else if (z>=am && z<a) {
            fct = (f[2].im+f[3].im*sqrt(a-z)+f[4].im*(a-z))/(z-x);
        }
        else if (z>a && z<=ap) {
            fct = (g[2].im+g[3].im*sqrt(z-a)+g[4].im*(z-a))/(z-x);
        }
        else {
            fct = (gsl_spline_eval(M_hat[1].im,z,acc2)-g[0].im-g[1].im*(z-a))/(pow(z-a,1.5)*(z-x));
        }
        return fct;
    }
    
    F_re_sing.function = &f_re_sing;
    F_im_sing.function = &f_im_sing;
    
    pts[0] = s0;
    pts[1] = a;
    pts[2] = L2;
    
    for (i=0; i<N; i++) {
        F_re_sing.params = &s[i];
        F_im_sing.params = &s[i];
        
        gsl_integration_qagp(&F_re_sing,pts,3,eps,eps,N_int,w,&result.re,&error);
        gsl_integration_qagp(&F_im_sing,pts,3,eps,eps,N_int,w,&result.im,&error);
        
        if (fabs(s[i])<1.e-5) {
            log_const = 1./L2+s[i]/(2.*L2*L2)+s[i]*s[i]/(3.*L2*L2*L2)+s[i]*s[i]*s[i]/(4.*L2*L2*L2*L2);
        }
        else {
            log_const = log(L2/(L2-s[i]))/s[i];
        }
        
        M_inhom[i].re = result.re+M_const.re*log_const;
        M_inhom[i].im = result.im+M_const.im*log_const;
        
        M_inhom[i].re += f[0].re*Q32[i].re+g[0].re*Q32[i].im;
        M_inhom[i].re += f[1].re*Q12[i].re+g[1].re*Q12[i].im;
        
        M_inhom[i].im += f[0].im*Q32[i].re+g[0].im*Q32[i].im;
        M_inhom[i].im += f[1].im*Q12[i].re+g[1].im*Q12[i].im;
        
        //printf("%.10e %.10e %.10e\n",s[i],M_inhom[i].re,M_inhom[i].im);
        
        M_inhom[i].re *= pow(s[i],(double)n)/M_PI;
        M_inhom[i].im *= pow(s[i],(double)n)/M_PI;
    }
    
    gsl_integration_workspace_free(w);
    gsl_interp_accel_free(acc1);
    gsl_interp_accel_free(acc2);
}

void inhomogenity_sqrt_cubed_complex(complex *M_inhom, complex_spline *M_hat, complex *f, complex *g, complex *s, complex *Q12_f, complex *Q32_f, complex *Q12_g, complex *Q32_g, double s0, double L2, int N, int n) {
    int i,N_int;
    double a,am,ap,eps,error,pts[3],r,phi;
    complex result,temp,M_const,log_const;
    
    a = METAP_M_MPION_SQUARED;
    
    N_int = 1500;
    eps = 1.0e-12;
    
    am = a-DELTA_A;
    ap = a+DELTA_A;
    
    gsl_interp_accel *acc1 = gsl_interp_accel_alloc();
    gsl_interp_accel *acc2 = gsl_interp_accel_alloc();
    
    gsl_integration_workspace *w = gsl_integration_workspace_alloc(N_int);
    
    gsl_function F_re_sing;
    gsl_function F_im_sing;
    
    M_const.re = gsl_spline_eval(M_hat[1].re,L2,acc2)*L2/pow(L2-a,1.5);
    M_const.im = gsl_spline_eval(M_hat[1].im,L2,acc2)*L2/pow(L2-a,1.5);
    
    double f_re_sing(double z, void * params) {
        complex x = *(complex *) params;
        double fct;
        if (z<am) {
            fct = ((gsl_spline_eval(M_hat[0].re,z,acc1)-f[0].re-f[1].re*(a-z))*(z-x.re)-(gsl_spline_eval(M_hat[0].im,z,acc1)-f[0].im-f[1].im*(a-z))*x.im)/(pow(a-z,1.5)*(pow(z-x.re,2.0)+x.im*x.im));
        }
        else if (z>=am && z<a) {
            fct = ((f[2].re+f[3].re*sqrt(a-z)+f[4].re*(a-z))*(z-x.re)-(f[2].im+f[3].im*sqrt(a-z)+f[4].im*(a-z))*x.im)/(pow(z-x.re,2.0)+x.im*x.im);
        }
        else if (z>a && z<=ap) {
            fct = ((g[2].re+g[3].re*sqrt(z-a)+g[4].re*(z-a))*(z-x.re)-(g[2].im+g[3].im*sqrt(z-a)+g[4].im*(z-a))*x.im)/(pow(z-x.re,2.0)+x.im*x.im);
        }
        else {
            fct = ((gsl_spline_eval(M_hat[1].re,z,acc2)-g[0].re-g[1].re*(z-a))*(z-x.re)-(gsl_spline_eval(M_hat[1].im,z,acc2)-g[0].im-g[1].im*(z-a))*x.im)/(pow(z-a,1.5)*(pow(z-x.re,2.0)+x.im*x.im));
        }
        return fct;
    }
    
    double f_im_sing(double z, void * params) {
        complex x = *(complex *) params;
        double fct;
        if (z<am) {
            fct = ((gsl_spline_eval(M_hat[0].im,z,acc1)-f[0].im-f[1].im*(a-z))*(z-x.re)+(gsl_spline_eval(M_hat[0].re,z,acc1)-f[0].re-f[1].re*(a-z))*x.im)/(pow(a-z,1.5)*(pow(z-x.re,2.0)+x.im*x.im));
        }
        else if (z>=am && z<a) {
            fct = ((f[2].im+f[3].im*sqrt(a-z)+f[4].im*(a-z))*(z-x.re)+(f[2].re+f[3].re*sqrt(a-z)+f[4].re*(a-z))*x.im)/(pow(z-x.re,2.0)+x.im*x.im);
        }
        else if (z>a && z<=ap) {
            fct = ((g[2].im+g[3].im*sqrt(z-a)+g[4].im*(z-a))*(z-x.re)+(g[2].re+g[3].re*sqrt(z-a)+g[4].re*(z-a))*x.im)/(pow(z-x.re,2.0)+x.im*x.im);
        }
        else {
            fct = ((gsl_spline_eval(M_hat[1].im,z,acc2)-g[0].im-g[1].im*(z-a))*(z-x.re)+(gsl_spline_eval(M_hat[1].re,z,acc2)-g[0].re-g[1].re*(z-a))*x.im)/(pow(z-a,1.5)*(pow(z-x.re,2.0)+x.im*x.im));
        }
        return fct;
    }
    
    F_re_sing.function = &f_re_sing;
    F_im_sing.function = &f_im_sing;
    
    pts[0] = s0;
    pts[1] = a;
    pts[2] = L2;
    
    for (i=0; i<N; i++) {
        F_re_sing.params = &s[i];
        F_im_sing.params = &s[i];
        
        gsl_integration_qagp(&F_re_sing,pts,3,eps,eps,N_int,w,&result.re,&error);
        gsl_integration_qagp(&F_im_sing,pts,3,eps,eps,N_int,w,&result.im,&error);
        
        log_const.re = 0.5*(M_PI*fabs(s[i].im)+2.*s[i].im*atan((s[i].re-L2)/s[i].im)-s[i].re*(log(s[i].im*s[i].im+pow(L2-s[i].re,2.))-2.*log(L2)))/(s[i].re*s[i].re+s[i].im*s[i].im);
        
        log_const.im = 0.5*(M_PI*s[i].re*sign(s[i].im)+2.*s[i].re*atan((s[i].re-L2)/s[i].im)+s[i].im*(log(s[i].im*s[i].im+pow(L2-s[i].re,2.))-2.*log(L2)))/(s[i].re*s[i].re+s[i].im*s[i].im);
        
        M_inhom[i].re = result.re+log_const.re*M_const.re-log_const.im*M_const.im;
        M_inhom[i].im = result.im+log_const.re*M_const.im+log_const.im*M_const.re;
        
        M_inhom[i].re += f[0].re*Q32_f[i].re-f[0].im*Q32_f[i].im+f[1].re*Q12_f[i].re-f[1].im*Q12_f[i].im+g[0].re*Q32_g[i].re-g[0].im*Q32_g[i].im+g[1].re*Q12_g[i].re-g[1].im*Q12_g[i].im;
        M_inhom[i].im += f[0].im*Q32_f[i].re+f[0].re*Q32_f[i].im+f[1].im*Q12_f[i].re+f[1].re*Q12_f[i].im+g[0].im*Q32_g[i].re+g[0].re*Q32_g[i].im+g[1].im*Q12_g[i].re+g[1].re*Q12_g[i].im;
        
        r = pow(sqrt(s[i].re*s[i].re+s[i].im*s[i].im),(double)n)/M_PI;
        phi = n*atan2(s[i].im,s[i].re);
        
        temp.re = r*(M_inhom[i].re*cos(phi)-M_inhom[i].im*sin(phi));
        temp.im = r*(M_inhom[i].re*sin(phi)+M_inhom[i].im*cos(phi));
        
        M_inhom[i].re = temp.re;
        M_inhom[i].im = temp.im;
        
        //printf("%d %.10e %.10e\n",i,M_inhom[i].re,M_inhom[i].im);
    }
    
    gsl_integration_workspace_free(w);
    gsl_interp_accel_free(acc1);
    gsl_interp_accel_free(acc2);
}

void build_inhomogenity_sqrt_cubed(complex **M_inhom, complex_spline *M_hat, complex *F, complex *G, double s0, double L2, double *s_cv_plus, double *s_cv_minus, double *s_below_cut, complex *s_complex, complex *R32_cv_plus, complex *Q12_cv_plus, complex *Q32_cv_plus, complex *R32_cv_minus, complex *Q12_cv_minus, complex *Q32_cv_minus, complex *Q12_below_cut, complex *Q32_below_cut, complex *Q12_complex_f, complex *Q32_complex_f, complex *Q12_complex_g, complex *Q32_complex_g, int *N, int n) {
    inhomogenity_sqrt_cubed_below_cut(M_inhom[3],M_hat,F,G,s_below_cut,Q12_below_cut,Q32_below_cut,s0,L2,N[3],n);
    inhomogenity_sqrt_cubed_cv_plus(M_inhom[0],M_hat,F,G,s_cv_plus,R32_cv_plus,Q12_cv_plus,Q32_cv_plus,s0,L2,1.,N[0],n);
    inhomogenity_sqrt_cubed_cv_minus(M_inhom[1],M_hat,F,G,s_cv_minus,R32_cv_minus,Q12_cv_minus,Q32_cv_minus,s0,L2,N[1],n);
    inhomogenity_sqrt_cubed_complex(M_inhom[2],M_hat,F,G,s_complex,Q12_complex_f,Q32_complex_f,Q12_complex_g,Q32_complex_g,s0,L2,N[2],n);
}

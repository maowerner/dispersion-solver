#include "Basic.h"
#include "InhomogenitySqrt.h"

//
void inhomogenity_sqrt_cv_plus(complex *M_inhom, complex_spline *M_hat, complex *f, complex *g, double *s, double s0, double a, double L2, double cut, int N, int n) {
    int i,N_int,Nm,Np;
    double am,ap,p,eps,error,interpol_re[4],interpol_im[4],s3,s2,pts_cv[2][3],pts_a[2][3],log_const;
    complex result,M,M_const,R12,Q12;
    
    N_int = 1500;
    eps = 1.0e-10;
    
    am = a-DELTA_A;
    ap = a+DELTA_A;
    
    gsl_interp_accel *acc1_re = gsl_interp_accel_alloc();
    gsl_interp_accel *acc1_im = gsl_interp_accel_alloc();
    gsl_interp_accel *acc2_re = gsl_interp_accel_alloc();
    gsl_interp_accel *acc2_im = gsl_interp_accel_alloc();
    
    gsl_integration_workspace *w = gsl_integration_workspace_alloc(N_int);
    
    gsl_function F_re_cv;
    gsl_function F_im_cv;
    gsl_function F_re_sing;
    gsl_function F_im_sing;
    
    M_const.re = gsl_spline_eval(M_hat[1].re,L2,acc2_re)*L2/sqrt(L2-a);
    M_const.im = gsl_spline_eval(M_hat[1].im,L2,acc2_im)*L2/sqrt(L2-a);
    
    double f_re_cv(double z, void * params) {
        double x = *(double *) params;
        double fct;
        if (z<a) {
            fct = (gsl_spline_eval(M_hat[0].re,z,acc1_re)-gsl_spline_eval(M_hat[0].re,x,acc1_re))/(sqrt(a-z)*(z-x));
        }
        else {
            fct = (gsl_spline_eval(M_hat[1].re,z,acc2_re)-gsl_spline_eval(M_hat[1].re,x,acc2_re))/(sqrt(z-a)*(z-x));
        }
        
        if (isnan(fct)==1) {
            printf("FATAL ERROR: In function inhomogenity_sqrt_cv_plus! fct=nan");
            exit(1);
        }
        
        return fct;
    }
    
    double f_im_cv(double z, void * params) {
        double x = *(double *) params;
        double fct;
        if (z<a) {
            fct = (gsl_spline_eval(M_hat[0].im,z,acc1_im)-gsl_spline_eval(M_hat[0].im,x,acc1_im))/(sqrt(a-z)*(z-x));
        }
        else {
            fct = (gsl_spline_eval(M_hat[1].im,z,acc2_im)-gsl_spline_eval(M_hat[1].im,x,acc2_im))/(sqrt(z-a)*(z-x));
        }
        return fct;
    }
    
    double f_re_sing(double z, void * params) {
        double x = *(double *) params;
        double fct;
        if (z<am) {
            fct = (gsl_spline_eval(M_hat[0].re,z,acc1_re)-f[0].re)/(sqrt(a-z)*(z-x));
        }
        else if (z>=am && z<a) {
            fct = (f[1].re+f[2].re*sqrt(a-z)+f[3].re*(a-z))/(z-x);
        }
        else if (z>a && z<=ap) {
            fct = (g[1].re+g[2].re*sqrt(z-a)+g[3].re*(z-a))/(z-x);
        }
        else {
            fct = (gsl_spline_eval(M_hat[1].re,z,acc2_re)-g[0].re)/(sqrt(z-a)*(z-x));
        }
        return fct;
    }
    
    double f_im_sing(double z, void * params) {
        double x = *(double *) params;
        double fct;
        if (z<am) {
            fct = (gsl_spline_eval(M_hat[0].im,z,acc1_im)-f[0].im)/(sqrt(a-z)*(z-x));
        }
        else if (z>=am && z<a) {
            fct = (f[1].im+f[2].im*sqrt(a-z)+f[3].im*(a-z))/(z-x);
        }
        else if (z>a && z<=ap) {
            fct = (g[1].im+g[2].im*sqrt(z-a)+g[3].im*(z-a))/(z-x);
        }
        else {
            fct = (gsl_spline_eval(M_hat[1].im,z,acc2_im)-g[0].im)/(sqrt(z-a)*(z-x));
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
        
        R12 = R_sqrt(s[i],s0,a,L2);
        Q12 = Q_sqrt(s[i],s0,a,L2);
        
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
            
            M.re = gsl_spline_eval(M_hat[0].re,s[i],acc1_re);
            M.im = gsl_spline_eval(M_hat[0].im,s[i],acc1_im);
            
            log_const = log(L2/(L2-s[i]))/s[i];
            
            M_inhom[i].re += M.re*R12.re-M.im*R12.im+f[0].re*Q12.re+g[0].re*Q12.im+M_const.re*log_const;
            M_inhom[i].im += M.re*R12.im+M.im*R12.re+f[0].im*Q12.re+g[0].im*Q12.im+M_const.im*log_const;
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
            
            M.re = gsl_spline_eval(M_hat[1].re,s[i],acc2_re);
            M.im = gsl_spline_eval(M_hat[1].im,s[i],acc2_im);
            
            log_const = log(L2/(L2-s[i]))/s[i];
            
            M_inhom[i].re += M.re*R12.re-M.im*R12.im+g[0].re*Q12.re+f[0].re*Q12.im+M_const.re*log_const;
            M_inhom[i].im += M.re*R12.im+M.im*R12.re+g[0].im*Q12.re+f[0].im*Q12.im+M_const.im*log_const;
        }
        
        else {
            Np = i+1;
        }
        
        //printf("%.10e %.10e %.10e\n",s[i],M_inhom[i].re,M_inhom[i].im);
        
        M_inhom[i].re *= pow(s[i],(double)n)/M_PI;
        M_inhom[i].im *= pow(s[i],(double)n)/M_PI;
    }
    
    cubic_spline_f01_df01(s[Nm],s[Np],M_inhom[Nm].re,M_inhom[Np].re,(M_inhom[Nm].re-M_inhom[Nm-1].re)/(s[Nm]-s[Nm-1]),(M_inhom[Np+1].re-M_inhom[Np].re)/(s[Np+1]-s[Np]),interpol_re);
    cubic_spline_f01_df01(s[Nm],s[Np],M_inhom[Nm].im,M_inhom[Np].im,(M_inhom[Nm].im-M_inhom[Nm-1].im)/(s[Nm]-s[Nm-1]),(M_inhom[Np+1].im-M_inhom[Np].im)/(s[Np+1]-s[Np]),interpol_im);
    
    for (i=Nm+1; i<Np; i++) {
        s2 = s[i]*s[i];
        s3 = s[i]*s2;
        M_inhom[i].re = interpol_re[0]*s3+interpol_re[1]*s2+interpol_re[2]*s[i]+interpol_re[3];
        M_inhom[i].im = interpol_im[0]*s3+interpol_im[1]*s2+interpol_im[2]*s[i]+interpol_im[3];
    }
    
    gsl_integration_workspace_free(w);
    gsl_interp_accel_free(acc1_re);
    gsl_interp_accel_free(acc1_im);
    gsl_interp_accel_free(acc2_re);
    gsl_interp_accel_free(acc2_im);
}

void inhomogenity_sqrt_cv_minus(complex *M_inhom, complex_spline *M_hat, complex *f, complex *g, double *s, double s0, double a, double L2, int N, int n) {
    int i,N_int;
    double am,ap,p,eps,error,pts_cv[3],pts_a[3],log_const;
    complex result,M,M_const,R12,Q12;
    
    N_int = 1500;
    eps = 1.0e-10;
    
    am = a-DELTA_A;
    ap = a+DELTA_A;
    
    gsl_interp_accel *acc1_re = gsl_interp_accel_alloc();
    gsl_interp_accel *acc1_im = gsl_interp_accel_alloc();
    gsl_interp_accel *acc2_re = gsl_interp_accel_alloc();
    gsl_interp_accel *acc2_im = gsl_interp_accel_alloc();
    
    gsl_integration_workspace *w = gsl_integration_workspace_alloc(N_int);
    
    M_const.re = gsl_spline_eval(M_hat[1].re,L2,acc2_re)*L2/sqrt(L2-a);
    M_const.im = gsl_spline_eval(M_hat[1].im,L2,acc2_im)*L2/sqrt(L2-a);
    
    gsl_function F_re_cv;
    gsl_function F_im_cv;
    gsl_function F_re_sing;
    gsl_function F_im_sing;
    
    double f_re_cv(double z, void * params) {
        double x = *(double *) params;
        double fct;
        if (z<a) {
            fct = (gsl_spline_eval(M_hat[0].re,z,acc1_re)-gsl_spline_eval(M_hat[0].re,x,acc1_re))/(sqrt(a-z)*(z-x));
        }
        else {
            fct = (gsl_spline_eval(M_hat[1].re,z,acc2_re)-gsl_spline_eval(M_hat[1].re,x,acc2_re))/(sqrt(z-a)*(z-x));
        }
        if (isnan(fct)==1) {
            //printf("re_cv %.10f %.10f %d\n",x,z,n);
        }
        return fct;
    }
    
    double f_im_cv(double z, void * params) {
        double x = *(double *) params;
        double fct;
        if (z<a) {
            fct = (gsl_spline_eval(M_hat[0].im,z,acc1_im)-gsl_spline_eval(M_hat[0].im,x,acc1_im))/(sqrt(a-z)*(z-x));
        }
        else {
            fct = (gsl_spline_eval(M_hat[1].im,z,acc2_im)-gsl_spline_eval(M_hat[1].im,x,acc2_im))/(sqrt(z-a)*(z-x));
        }
        if (isnan(fct)==1) {
            //printf("im_cv %.10f %.10f %d\n",x,z,n);
        }
        return fct;
    }
    
    double f_re_sing(double z, void * params) {
        double x = *(double *) params;
        double fct;
        if (z<am) {
            fct = (gsl_spline_eval(M_hat[0].re,z,acc1_re)-f[0].re)/(sqrt(a-z)*(z-x));
        }
        else if (z>=am && z<a) {
            fct = (f[1].re+f[2].re*sqrt(a-z)+f[3].re*(a-z))/(z-x);
        }
        else if (z>a && z<=ap) {
            fct = (g[1].re+g[2].re*sqrt(z-a)+g[3].re*(z-a))/(z-x);
        }
        else {
            fct = (gsl_spline_eval(M_hat[1].re,z,acc2_re)-g[0].re)/(sqrt(z-a)*(z-x));
        }
        if (isnan(fct)==1) {
            //printf("re_sin %.10f %.10f\n",x,z);
        }
        return fct;
    }
    
    double f_im_sing(double z, void * params) {
        double x = *(double *) params;
        double fct;
        if (z<am) {
            fct = (gsl_spline_eval(M_hat[0].im,z,acc1_im)-f[0].im)/(sqrt(a-z)*(z-x));
        }
        else if (z>=am && z<a) {
            fct = (f[1].im+f[2].im*sqrt(a-z)+f[3].im*(a-z))/(z-x);
        }
        else if (z>a && z<=ap) {
            fct = (g[1].im+g[2].im*sqrt(z-a)+g[3].im*(z-a))/(z-x);
        }
        else {
            fct = (gsl_spline_eval(M_hat[1].im,z,acc2_im)-g[0].im)/(sqrt(z-a)*(z-x));
        }
        if (isnan(fct)==1) {
            //printf("im_sin %.10f %.10f\n",x,z);
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
        
        F_re_cv.params = &s[i];
        F_im_cv.params = &s[i];
        
        F_re_sing.params = &s[i];
        F_im_sing.params = &s[i];
        
        pts_cv[1] = s[i];
        pts_cv[2] = p;
        
        pts_a[0] = p;
        
        gsl_integration_qagp(&F_re_cv,pts_cv,3,eps,eps,N_int,w,&result.re,&error);
        gsl_integration_qagp(&F_im_cv,pts_cv,3,eps,eps,N_int,w,&result.im,&error);
        
        M_inhom[i].re = result.re;
        M_inhom[i].im = result.im;
        
        gsl_integration_qagp(&F_re_sing,pts_a,3,eps,eps,N_int,w,&result.re,&error);
        gsl_integration_qagp(&F_im_sing,pts_a,3,eps,eps,N_int,w,&result.im,&error);
        
        M_inhom[i].re += result.re;
        M_inhom[i].im += result.im;
        
        M.re = gsl_spline_eval(M_hat[0].re,s[i],acc1_re);
        M.im = gsl_spline_eval(M_hat[0].im,s[i],acc1_im);
        
        log_const = log(L2/(L2-s[i]))/s[i];
        
        R12 = R_sqrt(s[i],s0,a,L2);
        Q12 = Q_sqrt(s[i],s0,a,L2);
        
        M_inhom[i].re += M.re*R12.re+M.im*R12.im+f[0].re*Q12.re+g[0].re*Q12.im+M_const.re*log_const;
        M_inhom[i].im += -M.re*R12.im+M.im*R12.re+f[0].im*Q12.re+g[0].im*Q12.im+M_const.im*log_const;
        
        M_inhom[i].re *= pow(s[i],(double)n)/M_PI;
        M_inhom[i].im *= pow(s[i],(double)n)/M_PI;
        
        if (isnan(M_inhom[i].re)==1) {
            printf("FATAL ERROR: In function inhomogenity_sqrt_cv_minus: M_inhom_%d[i].re=nan, step=%d\n",n,i);
            exit(1);
        }
        else if (isnan(M_inhom[i].im)==1) {
            printf("FATAL ERROR: In function inhomogenity_sqrt_cv_minus: M_inhom_%d[i].im=nan, step=%d\n",n,i);
            exit(1);
        }
        
        //printf("%.10e %.10e %.10e\n",s[i],M_inhom[i].re,M_inhom[i].im);
    }
    
    gsl_integration_workspace_free(w);
    gsl_interp_accel_free(acc1_re);
    gsl_interp_accel_free(acc1_im);
    gsl_interp_accel_free(acc2_re);
    gsl_interp_accel_free(acc2_im);
}

void inhomogenity_sqrt_below_cut(complex *M_inhom, complex_spline *M_hat, complex *f, complex *g, double *s, double s0, double a, double L2, int N, int n) {
    int i,N_int;
    double am,ap,eps,error,pts[3],log_const;
    complex result,M_const,Q12;
    
    double *M_re = (double *)malloc(N*sizeof(double));
    double *M_im = (double *)malloc(N*sizeof(double));
    
    N_int = 1500;
    eps = 1.0e-10;
    
    am = a-DELTA_A;
    ap = a+DELTA_A;
    
    gsl_interp_accel *acc1_re = gsl_interp_accel_alloc();
    gsl_interp_accel *acc1_im = gsl_interp_accel_alloc();
    gsl_interp_accel *acc2_re = gsl_interp_accel_alloc();
    gsl_interp_accel *acc2_im = gsl_interp_accel_alloc();
    
    gsl_integration_workspace *w = gsl_integration_workspace_alloc(N_int);
    
    M_const.re = gsl_spline_eval(M_hat[1].re,L2,acc2_re)*L2/sqrt(L2-a);
    M_const.im = gsl_spline_eval(M_hat[1].im,L2,acc2_im)*L2/sqrt(L2-a);
    
    gsl_function F_re_sing;
    gsl_function F_im_sing;
    
    double f_re_sing(double z, void * params) {
        double x = *(double *) params;
        double fct;
        if (z<am) {
            fct = (gsl_spline_eval(M_hat[0].re,z,acc1_re)-f[0].re)/(sqrt(a-z)*(z-x));
        }
        else if (z>=am && z<a) {
            fct = (f[1].re+f[2].re*sqrt(a-z)+f[3].re*(a-z))/(z-x);
        }
        else if (z>a && z<=ap) {
            fct = (g[1].re+g[2].re*sqrt(z-a)+g[3].re*(z-a))/(z-x);
        }
        else {
            fct = (gsl_spline_eval(M_hat[1].re,z,acc2_re)-g[0].re)/(sqrt(z-a)*(z-x));
        }
        
        if (isnan(fct)==1) {
            printf("%.3e %.3e\n",z,x);
        }
        
        return fct;
    }
    
    double f_im_sing(double z, void * params) {
        double x = *(double *) params;
        double fct;
        if (z<am) {
            fct = (gsl_spline_eval(M_hat[0].im,z,acc1_im)-f[0].im)/(sqrt(a-z)*(z-x));
        }
        else if (z>=am && z<a) {
            fct = (f[1].im+f[2].im*sqrt(a-z)+f[3].im*(a-z))/(z-x);
        }
        else if (z>a && z<=ap) {
            fct = (g[1].im+g[2].im*sqrt(z-a)+g[3].im*(z-a))/(z-x);
        }
        else {
            fct = (gsl_spline_eval(M_hat[1].im,z,acc2_im)-g[0].im)/(sqrt(z-a)*(z-x));
        }
        
        if (isnan(fct)==1) {
            printf("%.3e %.3e\n",z,x);
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
        
        M_inhom[i].re = result.re;
        M_inhom[i].im = result.im;
        
        if (fabs(s[i])<1.e-5) {
            log_const = 1./L2+s[i]/(2.*L2*L2)+s[i]*s[i]/(3.*L2*L2*L2)+s[i]*s[i]*s[i]/(4.*L2*L2*L2*L2);
        }
        else {
            log_const = log(L2/(L2-s[i]))/s[i];
        }
        
        Q12 = Q_sqrt(s[i],s0,a,L2);
        
        M_inhom[i].re += f[0].re*Q12.re+g[0].re*Q12.im+M_const.re*log_const;
        M_inhom[i].im += f[0].im*Q12.re+g[0].im*Q12.im+M_const.im*log_const;
        
        //printf("%.10e %.10e %.10e\n",s[i],M_inhom[i].re,M_inhom[i].im);
        
        M_inhom[i].re *= pow(s[i],(double)n)/M_PI;
        M_inhom[i].im *= pow(s[i],(double)n)/M_PI;
        
        if (isnan(M_inhom[i].re)==1) {
            printf("FATAL ERROR: In function inhomogenity_sqrt_below_cut: M_inhom_%d[i].re=nan, step=%d, s=%.10e\n",n,i,s[i]);
            exit(1);
        }
        else if (isnan(M_inhom[i].im)==1) {
            printf("FATAL ERROR: In function inhomogenity_sqrt_below_cut: M_inhom_%d[i].im=nan, step=%d, s=%10e\n",n,i,s[i]);
            exit(1);
        }
    }
    
    gsl_integration_workspace_free(w);
    gsl_interp_accel_free(acc1_re);
    gsl_interp_accel_free(acc1_im);
    gsl_interp_accel_free(acc2_re);
    gsl_interp_accel_free(acc2_im);
}

complex inhomogenity_sqrt_below_cut_single_value(complex_spline *M_hat, complex *f, complex *g, double s, double s0, double a, double L2, int n) {
    
    complex M_inhom;
    
    int i,N_int;
    double am,ap,eps,error,pts[3],log_const;
    complex result,M_const,Q12;
    
    //double *M_re = (double *)malloc(N*sizeof(double));
    //double *M_im = (double *)malloc(N*sizeof(double));
    
    N_int = 1500;
    eps = 1.0e-10;
    
    am = a-DELTA_A;
    ap = a+DELTA_A;
    
    gsl_interp_accel *acc1_re = gsl_interp_accel_alloc();
    gsl_interp_accel *acc1_im = gsl_interp_accel_alloc();
    gsl_interp_accel *acc2_re = gsl_interp_accel_alloc();
    gsl_interp_accel *acc2_im = gsl_interp_accel_alloc();
    
    gsl_integration_workspace *w = gsl_integration_workspace_alloc(N_int);
    
    M_const.re = gsl_spline_eval(M_hat[1].re,L2,acc2_re)*L2/sqrt(L2-a);
    M_const.im = gsl_spline_eval(M_hat[1].im,L2,acc2_im)*L2/sqrt(L2-a);
    
    gsl_function F_re_sing;
    gsl_function F_im_sing;
    
    double f_re_sing(double z, void * params) {
        double x = *(double *) params;
        double fct;
        if (z<am) {
            fct = (gsl_spline_eval(M_hat[0].re,z,acc1_re)-f[0].re)/(sqrt(a-z)*(z-x));
        }
        else if (z>=am && z<a) {
            fct = (f[1].re+f[2].re*sqrt(a-z)+f[3].re*(a-z))/(z-x);
        }
        else if (z>a && z<=ap) {
            fct = (g[1].re+g[2].re*sqrt(z-a)+g[3].re*(z-a))/(z-x);
        }
        else {
            fct = (gsl_spline_eval(M_hat[1].re,z,acc2_re)-g[0].re)/(sqrt(z-a)*(z-x));
        }
        
        if (isnan(fct)==1) {
            printf("%.3e %.3e\n",z,x);
        }
        
        return fct;
    }
    
    double f_im_sing(double z, void * params) {
        double x = *(double *) params;
        double fct;
        if (z<am) {
            fct = (gsl_spline_eval(M_hat[0].im,z,acc1_im)-f[0].im)/(sqrt(a-z)*(z-x));
        }
        else if (z>=am && z<a) {
            fct = (f[1].im+f[2].im*sqrt(a-z)+f[3].im*(a-z))/(z-x);
        }
        else if (z>a && z<=ap) {
            fct = (g[1].im+g[2].im*sqrt(z-a)+g[3].im*(z-a))/(z-x);
        }
        else {
            fct = (gsl_spline_eval(M_hat[1].im,z,acc2_im)-g[0].im)/(sqrt(z-a)*(z-x));
        }
        
        if (isnan(fct)==1) {
            printf("%.3e %.3e\n",z,x);
        }
        
        return fct;
    }
    
    F_re_sing.function = &f_re_sing;
    F_im_sing.function = &f_im_sing;
    
    pts[0] = s0;
    pts[1] = a;
    pts[2] = L2;
    
    F_re_sing.params = &s;
    F_im_sing.params = &s;
    
    gsl_integration_qagp(&F_re_sing,pts,3,eps,eps,N_int,w,&result.re,&error);
    gsl_integration_qagp(&F_im_sing,pts,3,eps,eps,N_int,w,&result.im,&error);
    
    M_inhom.re = result.re;
    M_inhom.im = result.im;
    
    if (fabs(s)<1.e-5) {
        log_const = 1./L2+s/(2.*L2*L2)+s*s/(3.*L2*L2*L2)+s*s*s/(4.*L2*L2*L2*L2);
    }
    else {
        log_const = log(L2/(L2-s))/s;
    }
    
    Q12 = Q_sqrt(s,s0,a,L2);
    
    M_inhom.re += f[0].re*Q12.re+g[0].re*Q12.im+M_const.re*log_const;
    M_inhom.im += f[0].im*Q12.re+g[0].im*Q12.im+M_const.im*log_const;
    
    //printf("%.10e %.10e %.10e\n",s[i],M_inhom[i].re,M_inhom[i].im);
    
    M_inhom.re *= 1./M_PI;
    M_inhom.im *= 1./M_PI;
    
    if (isnan(M_inhom.re)==1) {
        printf("FATAL ERROR: In function inhomogenity_sqrt_below_cut: M_inhom_%d[i].re=nan, step=%d, s=%.10e\n",n,i,s);
        exit(1);
    }
    else if (isnan(M_inhom.im)==1) {
        printf("FATAL ERROR: In function inhomogenity_sqrt_below_cut: M_inhom_%d[i].im=nan, step=%d, s=%10e\n",n,i,s);
        exit(1);
    }
    
    gsl_integration_workspace_free(w);
    gsl_interp_accel_free(acc1_re);
    gsl_interp_accel_free(acc1_im);
    gsl_interp_accel_free(acc2_re);
    gsl_interp_accel_free(acc2_im);
    
    return M_inhom;
}

void inhomogenity_sqrt_complex(complex *M_inhom, complex_spline *M_hat, complex *f, complex *g, complex *s, char plusminus, double s0, double a, double L2, int N, int n) {
    int i,N_int;
    double am,ap,eps,error,pts[3],r,phi;
    complex result,temp,M_const,log_const,Q12_f,Q12_g;
    
    N_int = 1500;
    eps = 1.0e-10;
    
    am = a-DELTA_A;
    ap = a+DELTA_A;
    
    gsl_interp_accel *acc1 = gsl_interp_accel_alloc();
    gsl_interp_accel *acc2 = gsl_interp_accel_alloc();
    
    gsl_integration_workspace *w = gsl_integration_workspace_alloc(N_int);
    
    M_const.re = gsl_spline_eval(M_hat[1].re,L2,acc2)*L2/sqrt(L2-a);
    M_const.im = gsl_spline_eval(M_hat[1].im,L2,acc2)*L2/sqrt(L2-a);
    
    gsl_function F_re_sing;
    gsl_function F_im_sing;
    
    double f_re_sing(double z, void * params) {
        complex x = *(complex *) params;
        double fct;
        if (z<am) {
            fct = ((gsl_spline_eval(M_hat[0].re,z,acc1)-f[0].re)*(z-x.re)-(gsl_spline_eval(M_hat[0].im,z,acc1)-f[0].im)*x.im)/(sqrt(a-z)*(pow(z-x.re,2.0)+x.im*x.im));
        }
        else if (z>=am && z<a) {
            fct = ((f[1].re+f[2].re*sqrt(a-z)+f[3].re*(a-z))*(z-x.re)-(f[1].im+f[2].im*sqrt(a-z)+f[3].im*(a-z))*x.im)/(pow(z-x.re,2.0)+x.im*x.im);
        }
        else if (z>a && z<=ap) {
            fct = ((g[1].re+g[2].re*sqrt(z-a)+g[3].re*(z-a))*(z-x.re)-(g[1].im+g[2].im*sqrt(z-a)+g[3].im*(z-a))*x.im)/(pow(z-x.re,2.0)+x.im*x.im);
        }
        else {
            fct = ((gsl_spline_eval(M_hat[1].re,z,acc2)-g[0].re)*(z-x.re)-(gsl_spline_eval(M_hat[1].im,z,acc2)-g[0].im)*x.im)/(sqrt(z-a)*(pow(z-x.re,2.0)+x.im*x.im));
        }
        return fct;
    }
    
    double f_im_sing(double z, void * params) {
        complex x = *(complex *) params;
        double fct;
        if (z<am) {
            fct = ((gsl_spline_eval(M_hat[0].im,z,acc1)-f[0].im)*(z-x.re)+(gsl_spline_eval(M_hat[0].re,z,acc1)-f[0].re)*x.im)/(sqrt(a-z)*(pow(z-x.re,2.0)+x.im*x.im));
        }
        else if (z>=am && z<a) {
            fct = ((f[1].im+f[2].im*sqrt(a-z)+f[3].im*(a-z))*(z-x.re)+(f[1].re+f[2].re*sqrt(a-z)+f[3].re*(a-z))*x.im)/(pow(z-x.re,2.0)+x.im*x.im);
        }
        else if (z>a && z<=ap) {
            fct = ((g[1].im+g[2].im*sqrt(z-a)+g[3].im*(z-a))*(z-x.re)+(g[1].re+g[2].re*sqrt(z-a)+g[3].re*(z-a))*x.im)/(pow(z-x.re,2.0)+x.im*x.im);
        }
        else {
            fct = ((gsl_spline_eval(M_hat[1].im,z,acc2)-g[0].im)*(z-x.re)+(gsl_spline_eval(M_hat[1].re,z,acc2)-g[0].re)*x.im)/(sqrt(z-a)*(pow(z-x.re,2.0)+x.im*x.im));
        }
        return fct;
    }
    
    F_re_sing.function = &f_re_sing;
    F_im_sing.function = &f_im_sing;
    
    pts[0] = s0;
    pts[1] = a;
    pts[2] = L2;
    
    for (i=0; i<N; i++) {
        switch (plusminus) {
            case '+':
                break;
                
            case '-':
                s[i].im *= -1.;
                break;
                
            default:
                printf("FATAL ERROR: In function inhomogenity_sqrt_complex, plusminus has to be '+' or '-', '%c' is unknown!",plusminus);
                exit(1);
                break;
        }
        
        F_re_sing.params = &s[i];
        F_im_sing.params = &s[i];
        
        gsl_integration_qagp(&F_re_sing,pts,3,eps,eps,N_int,w,&result.re,&error);
        gsl_integration_qagp(&F_im_sing,pts,3,eps,eps,N_int,w,&result.im,&error);
        
        log_const.re = 0.5*(M_PI*fabs(s[i].im)+2.*s[i].im*atan((s[i].re-L2)/s[i].im)-s[i].re*(log(s[i].im*s[i].im+pow(L2-s[i].re,2.))-2.*log(L2)))/(s[i].re*s[i].re+s[i].im*s[i].im);
        
        log_const.im = 0.5*(M_PI*s[i].re*sign(s[i].im)+2.*s[i].re*atan((s[i].re-L2)/s[i].im)+s[i].im*(log(s[i].im*s[i].im+pow(L2-s[i].re,2.))-2.*log(L2)))/(s[i].re*s[i].re+s[i].im*s[i].im);
        
        M_inhom[i].re = result.re+log_const.re*M_const.re-log_const.im*M_const.im;
        M_inhom[i].im = result.im+log_const.re*M_const.im+log_const.im*M_const.re;
        
        Q12_f = Q_sqrt_complex_f(s[i],s0,a);
        Q12_g = Q_sqrt_complex_g(s[i],a,L2);
        
        M_inhom[i].re += f[0].re*Q12_f.re-f[0].im*Q12_f.im+g[0].re*Q12_g.re-g[0].im*Q12_g.im;
        M_inhom[i].im += f[0].im*Q12_f.re+f[0].re*Q12_f.im+g[0].im*Q12_g.re+g[0].re*Q12_g.im;
        
        r = pow(sqrt(s[i].re*s[i].re+s[i].im*s[i].im),(double)n)/M_PI;
        phi = n*atan2(s[i].im,s[i].re);
        
        temp.re = r*(M_inhom[i].re*cos(phi)-M_inhom[i].im*sin(phi));
        temp.im = r*(M_inhom[i].re*sin(phi)+M_inhom[i].im*cos(phi));
        
        M_inhom[i].re = temp.re;
        M_inhom[i].im = temp.im;
        
        if (isnan(M_inhom[i].re)==1) {
            printf("FATAL ERROR: In function inhomogenity_sqrt_complex: M_inhom_%d[i].re=nan, step=%d, s.re=%.10e s.im=%.10e\n",n,i,s[i].re,s[i].im);
            exit(1);
        }
        else if (isnan(M_inhom[i].im)==1) {
            printf("FATAL ERROR: In function inhomogenity_sqrt_complex: M_inhom_%d[i].im=nan, step=%d, s.re=%.10e s.im=%.10e\n",n,i,s[i].re,s[i].im);
            exit(1);
        }
        
        switch (plusminus) {
            case '+':
                break;
                
            case '-':
                s[i].im *= -1.;
                break;
                
            default:
                printf("FATAL ERROR: In function inhomogenity_sqrt_complex, plusminus has to be '+' or '-', '%c' is unknown!",plusminus);
                exit(1);
                break;
        }
        
        //printf("%d %.10e %.10e\n",i,M_inhom[i].re,M_inhom[i].im);
    }
    
    gsl_integration_workspace_free(w);
    gsl_interp_accel_free(acc1);
    gsl_interp_accel_free(acc2);
}

void build_inhomogenity_sqrt(complex **M_inhom, complex_spline *M_hat, complex *F, complex *G, double s0, double a, double L2, double *s_cv_plus, double *s_cv_minus, double *s_below_cut, complex **s_cmp, int *N, int n) {
    inhomogenity_sqrt_below_cut(M_inhom[2],M_hat,F,G,s_below_cut,s0,a,L2,N[2],n);
    inhomogenity_sqrt_cv_plus(M_inhom[0],M_hat,F,G,s_cv_plus,s0,a,L2,1.,N[0],n);
    inhomogenity_sqrt_cv_minus(M_inhom[1],M_hat,F,G,s_cv_minus,s0,a,L2,N[1],n);
    inhomogenity_sqrt_complex(M_inhom[3],M_hat,F,G,s_cmp[0],'-',s0,a,L2,N[3],n);
    inhomogenity_sqrt_complex(M_inhom[4],M_hat,F,G,s_cmp[1],'-',s0,a,L2,N[3],n);
    inhomogenity_sqrt_complex(M_inhom[5],M_hat,F,G,s_cmp[0],'+',s0,a,L2,N[3],n);
    inhomogenity_sqrt_complex(M_inhom[6],M_hat,F,G,s_cmp[1],'+',s0,a,L2,N[3],n);
}

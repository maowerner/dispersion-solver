#include "Basic.h"
#include "Singularity.h"
#include "InhomogenitySqrt.h"
#include "Etapi_disc.h"

void build_inhomogenity_integrand_etapi_disc(complex_spline M0_finite, complex_spline M2_finite, complex_spline *M0_sing, complex_spline *M2_sing, complex *cf0, complex *cg0, complex *cf2, complex *cg2, complex_spline f0etapi, complex_spline f2etapi, complex_spline Metapi, complex_spline *Metapi_tilde, gsl_spline *delta0, gsl_spline *delta2, double *abs_omnes0, double *abs_omnes2, double *s, double s0, double L2, int *N, int Nsin, int n0, int n2) {
    int i,j,n,m,N1,N2;
    double si,sf,s_am,s_ap,l,d0,d2,c0,c2,*M0_re,*M0_im,*M2_re,*M2_im,*M0_tilde_re,*M0_tilde_im,*M2_tilde_re,*M2_tilde_im,s_f[5],s_g[5];
    complex temp0,temp2,f0[4],g0[4],f2[4],g2[4];
    
    n = N[0]+N[1];
    
    M0_re = (double *)malloc(n*sizeof(double));
    M0_im = (double *)malloc(n*sizeof(double));
    
    M2_re = (double *)malloc(n*sizeof(double));
    M2_im = (double *)malloc(n*sizeof(double));
    
    M0_tilde_re = (double *)malloc(n*sizeof(double));
    M0_tilde_im = (double *)malloc(n*sizeof(double));
    
    M2_tilde_re = (double *)malloc(n*sizeof(double));
    M2_tilde_im = (double *)malloc(n*sizeof(double));
    
    gsl_interp_accel *acc_d0 = gsl_interp_accel_alloc();
    gsl_interp_accel *acc_d2 = gsl_interp_accel_alloc();
    gsl_interp_accel *acc_f0 = gsl_interp_accel_alloc();
    gsl_interp_accel *acc_f2 = gsl_interp_accel_alloc();
    gsl_interp_accel *acc_M = gsl_interp_accel_alloc();
    gsl_interp_accel *acc_M_tilde_low = gsl_interp_accel_alloc();
    gsl_interp_accel *acc_M_tilde_high = gsl_interp_accel_alloc();
    
    for (i=0; i<N[0]; i++) {
        l = lambda_etapi(s[i]);
        d0 = gsl_spline_eval(delta0,s[i],acc_d0);
        d2 = gsl_spline_eval(delta2,s[i],acc_d2);
        
        c0 = l/(48.*M_PI*abs_omnes0[i]*pow(s[i],n0-1.));
        c2 = l/(32.*M_PI*abs_omnes2[i]*pow(s[i],n2-1.));
        
        temp0.re = gsl_spline_eval(f0etapi.re,s[i],acc_f0)*gsl_spline_eval(Metapi.re,s[i],acc_M)+gsl_spline_eval(f0etapi.im,s[i],acc_f0)*gsl_spline_eval(Metapi.im,s[i],acc_M);
        temp0.im = gsl_spline_eval(f0etapi.re,s[i],acc_f0)*gsl_spline_eval(Metapi.im,s[i],acc_M)-gsl_spline_eval(f0etapi.im,s[i],acc_f0)*gsl_spline_eval(Metapi.re,s[i],acc_M);
        temp2.re = gsl_spline_eval(f2etapi.re,s[i],acc_f2)*gsl_spline_eval(Metapi.re,s[i],acc_M)+gsl_spline_eval(f2etapi.im,s[i],acc_f2)*gsl_spline_eval(Metapi.im,s[i],acc_M);
        temp2.im = gsl_spline_eval(f2etapi.re,s[i],acc_f2)*gsl_spline_eval(Metapi.im,s[i],acc_M)-gsl_spline_eval(f2etapi.im,s[i],acc_f2)*gsl_spline_eval(Metapi.re,s[i],acc_M);
        
        M0_re[i] = (temp0.re*cos(d0)-temp0.im*sin(d0))*c0;
        M0_im[i] = (temp0.re*sin(d0)+temp0.im*cos(d0))*c0;
        M2_re[i] = (temp2.re*cos(d0)-temp2.im*sin(d0))*c0;
        M2_im[i] = (temp2.re*sin(d0)+temp2.im*cos(d0))*c0;
        
        c0 /= s[i];
        c2 /= s[i];
        
        temp0.re = gsl_spline_eval(f0etapi.re,s[i],acc_f0)*gsl_spline_eval(Metapi_tilde[0].re,s[i],acc_M_tilde_low)+gsl_spline_eval(f0etapi.im,s[i],acc_f0)*gsl_spline_eval(Metapi_tilde[0].im,s[i],acc_M_tilde_low);
        temp0.im = gsl_spline_eval(f0etapi.re,s[i],acc_f0)*gsl_spline_eval(Metapi_tilde[0].im,s[i],acc_M_tilde_low)-gsl_spline_eval(f0etapi.im,s[i],acc_f0)*gsl_spline_eval(Metapi_tilde[0].re,s[i],acc_M_tilde_low);
        temp2.re = gsl_spline_eval(f2etapi.re,s[i],acc_f2)*gsl_spline_eval(Metapi_tilde[0].re,s[i],acc_M_tilde_low)+gsl_spline_eval(f2etapi.im,s[i],acc_f2)*gsl_spline_eval(Metapi_tilde[0].im,s[i],acc_M_tilde_low);
        temp2.im = gsl_spline_eval(f2etapi.re,s[i],acc_f2)*gsl_spline_eval(Metapi_tilde[0].im,s[i],acc_M_tilde_low)-gsl_spline_eval(f2etapi.im,s[i],acc_f2)*gsl_spline_eval(Metapi_tilde[0].re,s[i],acc_M_tilde_low);
        
        M0_tilde_re[i] = (temp0.re*cos(d0)-temp0.im*sin(d0))*c0;
        M0_tilde_im[i] = (temp0.re*sin(d0)+temp0.im*cos(d0))*c0;
        M2_tilde_re[i] = (temp2.re*cos(d0)-temp2.im*sin(d0))*c0;
        M2_tilde_im[i] = (temp2.re*sin(d0)+temp2.im*cos(d0))*c0;
    }
    for (i=N[0]; i<n; i++) {
        l = lambda_etapi(s[i]);
        d0 = gsl_spline_eval(delta0,s[i],acc_d0);
        d2 = gsl_spline_eval(delta2,s[i],acc_d2);
        
        c0 = l/(48.*M_PI*abs_omnes0[i]*pow(s[i],n0-1.));
        c2 = l/(32.*M_PI*abs_omnes2[i]*pow(s[i],n2-1.));
        
        temp0.re = gsl_spline_eval(f0etapi.re,s[i],acc_f0)*gsl_spline_eval(Metapi.re,s[i],acc_M)+gsl_spline_eval(f0etapi.im,s[i],acc_f0)*gsl_spline_eval(Metapi.im,s[i],acc_M);
        temp0.im = gsl_spline_eval(f0etapi.re,s[i],acc_f0)*gsl_spline_eval(Metapi.im,s[i],acc_M)-gsl_spline_eval(f0etapi.im,s[i],acc_f0)*gsl_spline_eval(Metapi.re,s[i],acc_M);
        temp2.re = gsl_spline_eval(f2etapi.re,s[i],acc_f2)*gsl_spline_eval(Metapi.re,s[i],acc_M)+gsl_spline_eval(f2etapi.im,s[i],acc_f2)*gsl_spline_eval(Metapi.im,s[i],acc_M);
        temp2.im = gsl_spline_eval(f2etapi.re,s[i],acc_f2)*gsl_spline_eval(Metapi.im,s[i],acc_M)-gsl_spline_eval(f2etapi.im,s[i],acc_f2)*gsl_spline_eval(Metapi.re,s[i],acc_M);
        
        M0_re[i] = (temp0.re*cos(d0)-temp0.im*sin(d0))*c0;
        M0_im[i] = (temp0.re*sin(d0)+temp0.im*cos(d0))*c0;
        M2_re[i] = (temp2.re*cos(d0)-temp2.im*sin(d0))*c0;
        M2_im[i] = (temp2.re*sin(d0)+temp2.im*cos(d0))*c0;
        
        c0 /= s[i];
        c2 /= s[i];
        
        temp0.re = gsl_spline_eval(f0etapi.re,s[i],acc_f0)*gsl_spline_eval(Metapi_tilde[1].re,s[i],acc_M_tilde_high)+gsl_spline_eval(f0etapi.im,s[i],acc_f0)*gsl_spline_eval(Metapi_tilde[1].im,s[i],acc_M_tilde_high);
        temp0.im = gsl_spline_eval(f0etapi.re,s[i],acc_f0)*gsl_spline_eval(Metapi_tilde[1].im,s[i],acc_M_tilde_high)-gsl_spline_eval(f0etapi.im,s[i],acc_f0)*gsl_spline_eval(Metapi_tilde[1].re,s[i],acc_M_tilde_high);
        temp2.re = gsl_spline_eval(f2etapi.re,s[i],acc_f2)*gsl_spline_eval(Metapi_tilde[1].re,s[i],acc_M_tilde_high)+gsl_spline_eval(f2etapi.im,s[i],acc_f2)*gsl_spline_eval(Metapi_tilde[1].im,s[i],acc_M_tilde_high);
        temp2.im = gsl_spline_eval(f2etapi.re,s[i],acc_f2)*gsl_spline_eval(Metapi_tilde[1].im,s[i],acc_M_tilde_high)-gsl_spline_eval(f2etapi.im,s[i],acc_f2)*gsl_spline_eval(Metapi_tilde[1].re,s[i],acc_M_tilde_high);
        
        M0_tilde_re[i] = (temp0.re*cos(d0)-temp0.im*sin(d0))*c0;
        M0_tilde_im[i] = (temp0.re*sin(d0)+temp0.im*cos(d0))*c0;
        M2_tilde_re[i] = (temp2.re*cos(d0)-temp2.im*sin(d0))*c0;
        M2_tilde_im[i] = (temp2.re*sin(d0)+temp2.im*cos(d0))*c0;
        
        //printf("%+.10e %+.10e %+.10e %+.10e %+.10e %+.10e %+.10e %+.10e\n",s[i],d0,abs_omnes0[i],gsl_spline_eval(f0etapi.re,s[i],acc_f0),gsl_spline_eval(f0etapi.im,s[i],acc_f0),gsl_spline_eval(Metapi.re,s[i],acc_M),gsl_spline_eval(Metapi.im,s[i],acc_M));
    }
    
    gsl_interp_accel_free(acc_d0);
    gsl_interp_accel_free(acc_d2);
    gsl_interp_accel_free(acc_f0);
    gsl_interp_accel_free(acc_f2);
    gsl_interp_accel_free(acc_M);
    gsl_interp_accel_free(acc_M_tilde_low);
    gsl_interp_accel_free(acc_M_tilde_high);
    
    m = 0;
    N1 = N[0]-1;
    N2 = N[0];
    for (i=0; i<4; i++) {
        f0[i].re = M0_tilde_re[N1-m];
        f0[i].im = M0_tilde_im[N1-m];
        
        f2[i].re = M2_tilde_re[N1-m];
        f2[i].im = M2_tilde_im[N1-m];
        
        g0[i].re = M0_tilde_re[N2+i];
        g0[i].im = M0_tilde_im[N2+i];
        
        g2[i].re = M2_tilde_re[N2+i];
        g2[i].im = M2_tilde_im[N2+i];
        
        s_f[i] = sqrt(fabs(s[N1-m]-METAP_M_MPION_SQUARED));
        s_g[i] = sqrt(fabs(s[N2+i]-METAP_M_MPION_SQUARED));
        m++;
        if (i==0) {
            N1 -= Nsin;
            N2 += Nsin;
        }
    }
    
    s_f[i] = sqrt(fabs(s[N1-m]-METAP_M_MPION_SQUARED));
    s_g[i] = sqrt(fabs(s[N2+i]-METAP_M_MPION_SQUARED));
    
    singularity_sqrt(f0,s_f,cf0);
    singularity_sqrt(g0,s_g,cg0);
    
    singularity_sqrt(f2,s_f,cf2);
    singularity_sqrt(g2,s_g,cg2);
    
    for (i=(N[0]-Nsin); i<N[0]; i++) {
        si = sqrt(METAP_M_MPION_SQUARED-s[i]);
        M0_tilde_re[i] = cf0[0].re+cf0[1].re*si+cf0[2].re*si*si+cf0[3].re*si*si*si;
        M0_tilde_im[i] = cf0[0].im+cf0[1].im*si+cf0[2].im*si*si+cf0[3].im*si*si*si;
        M2_tilde_re[i] = cf2[0].re+cf2[1].re*si+cf2[2].re*si*si+cf2[3].re*si*si*si;
        M2_tilde_im[i] = cf2[0].im+cf2[1].im*si+cf2[2].im*si*si+cf2[3].im*si*si*si;
    }
    for (i=N[0]; i<(N[0]+Nsin); i++) {
        si = sqrt(s[i]-METAP_M_MPION_SQUARED);
        M0_tilde_re[i] = cg0[0].re+cg0[1].re*si+cg0[2].re*si*si+cg0[3].re*si*si*si;
        M0_tilde_im[i] = cg0[0].im+cg0[1].im*si+cg0[2].im*si*si+cg0[3].im*si*si*si;
        M2_tilde_re[i] = cg2[0].re+cg2[1].re*si+cg2[2].re*si*si+cg2[3].re*si*si*si;
        M2_tilde_im[i] = cg2[0].im+cg2[1].im*si+cg2[2].im*si*si+cg2[3].im*si*si*si;
    }
    
//    printf("a=%+.10e\n",METAP_M_MPION_SQUARED);
//    printf("f0_re(s)=%+.10e%+.10e*sqrt(a-s)%+.10e*(a-s)%+.10e*sqrt(a-s)*(a-s)\n",cf0[0].re,cf0[1].re,cf0[2].re,cf0[3].re);
//    printf("f0_im(s)=%+.10e%+.10e*sqrt(a-s)%+.10e*(a-s)%+.10e*sqrt(a-s)*(a-s)\n",cf0[0].im,cf0[1].im,cf0[2].im,cf0[3].im);
//    printf("g0_re(s)=%+.10e%+.10e*sqrt(s-a)%+.10e*(s-a)%+.10e*sqrt(s-a)*(s-a)\n",cg0[0].re,cg0[1].re,cg0[2].re,cg0[3].re);
//    printf("g0_im(s)=%+.10e%+.10e*sqrt(s-a)%+.10e*(s-a)%+.10e*sqrt(s-a)*(s-a)\n",cg0[0].im,cg0[1].im,cg0[2].im,cg0[3].im);
//    printf("f2_re(s)=%+.10e%+.10e*sqrt(a-s)%+.10e*(a-s)%+.10e*sqrt(a-s)*(a-s)\n",cf2[0].re,cf2[1].re,cf2[2].re,cf2[3].re);
//    printf("f2_im(s)=%+.10e%+.10e*sqrt(a-s)%+.10e*(a-s)%+.10e*sqrt(a-s)*(a-s)\n",cf2[0].im,cf2[1].im,cf2[2].im,cf2[3].im);
//    printf("g2_re(s)=%+.10e%+.10e*sqrt(s-a)%+.10e*(s-a)%+.10e*sqrt(s-a)*(s-a)\n",cg2[0].re,cg2[1].re,cg2[2].re,cg2[3].re);
//    printf("g2_im(s)=%+.10e%+.10e*sqrt(s-a)%+.10e*(s-a)%+.10e*sqrt(s-a)*(s-a)\n",cg2[0].im,cg2[1].im,cg2[2].im,cg2[3].im);

    for (i=0; i<n; i++) {
        //printf("%+.10e %+.10e %+.10e %+.10e %+.10e %+.10e %+.10e %+.10e %+.10e\n",s[i],M0_re[i],M0_im[i],M2_re[i],M2_im[i],M0_tilde_re[i],M0_tilde_im[i],M2_tilde_re[i],M2_tilde_im[i]);
        if (isnan(M0_re[i])==1 || isnan(M0_im[i])==1) {
            printf("FATAL ERROR: M0_finite.re=%g, M0_finite.re=%g, step=%d\n",M0_re[i],M0_im[i],i);
            exit(1);
        }
        if (isnan(M2_re[i])==1 || isnan(M2_im[i])==1) {
            printf("FATAL ERROR: M2_finite.re=%g, M2_finite.re=%g, step=%d\n",M2_re[i],M2_im[i],i);
            exit(1);
        }
        if (isnan(M0_tilde_re[i])==1 || isnan(M0_tilde_im[i])==1) {
            printf("FATAL ERROR: M0_sing.re=%g, M0_sing.re=%g, step=%d\n",M0_tilde_re[i],M0_tilde_im[i],i);
            exit(1);
        }
        if (isnan(M2_tilde_re[i])==1 || isnan(M2_tilde_im[i])==1) {
            printf("FATAL ERROR: M2_sing.re=%g, M2_sing.re=%g, step=%d\n",M2_tilde_re[i],M2_tilde_im[i],i);
            exit(1);
        }
    }
    
    
    si = s[0];
    s_am = s[N[0]-1];
    s_ap = s[N[0]];
    sf = s[n-1];
    
    s[0] = s0;
    s[n-1] = L2;
    
    gsl_spline_init(M0_finite.re,s,M0_re,n);
    gsl_spline_init(M0_finite.im,s,M0_im,n);
    gsl_spline_init(M2_finite.re,s,M2_re,n);
    gsl_spline_init(M2_finite.im,s,M2_im,n);
    
    s[N[0]-1] = METAP_M_MPION_SQUARED;
    s[N[0]] = METAP_M_MPION_SQUARED;
    
    gsl_spline_init(M0_sing[0].re,s,M0_tilde_re,N[0]);
    gsl_spline_init(M0_sing[0].im,s,M0_tilde_im,N[0]);
    gsl_spline_init(M2_sing[0].re,s,M2_tilde_re,N[0]);
    gsl_spline_init(M2_sing[0].im,s,M2_tilde_im,N[0]);
    
    gsl_spline_init(M0_sing[1].re,s+N[0],M0_tilde_re+N[0],N[1]);
    gsl_spline_init(M0_sing[1].im,s+N[0],M0_tilde_im+N[0],N[1]);
    gsl_spline_init(M2_sing[1].re,s+N[0],M2_tilde_re+N[0],N[1]);
    gsl_spline_init(M2_sing[1].im,s+N[0],M2_tilde_im+N[0],N[1]);
    
    s[0] = si;
    s[N[0]-1] = s_am;
    s[N[0]] = s_ap;
    s[n-1] = sf;

    free(M0_re);
    free(M0_im);
    free(M2_re);
    free(M2_im);
    
    free(M0_tilde_re);
    free(M0_tilde_im);
    free(M2_tilde_re);
    free(M2_tilde_im);
}

//etapi_disc_finite function for given s in [s_etapi,L2] above the cut
void etapi_disc_finite_cv_plus(complex *M_etapi, complex_spline M_etapi_hat, double s0, double L2, double *s, int N, int n) {
    int i,N_int;
    double error,eps,pts[3],log_1,log_2;
    complex temp,M,Mconst;
    
    N_int  = 1500;
    eps = 1.0e-14;
    
    gsl_interp_accel *acc = gsl_interp_accel_alloc();
    gsl_integration_workspace *w = gsl_integration_workspace_alloc(N_int);
    gsl_function F_re;
    gsl_function F_im;
    
    double f_re (double z, void * params) {
        double x = *(double *) params;
        double fct = (gsl_spline_eval(M_etapi_hat.re,z,acc)-gsl_spline_eval(M_etapi_hat.re,x,acc))/(z*(z-x));
        return fct;
    }
    double f_im (double z, void * params) {
        double x = *(double *) params;
        double fct = (gsl_spline_eval(M_etapi_hat.im,z,acc)-gsl_spline_eval(M_etapi_hat.im,x,acc))/(z*(z-x));
        return fct;
    }
    
    F_re.function = &f_re;
    F_im.function = &f_im;
    
    Mconst.re = gsl_spline_eval(M_etapi_hat.re,L2,acc);
    Mconst.im = gsl_spline_eval(M_etapi_hat.im,L2,acc);
    
    pts[0] = s0;
    pts[2] = L2;
    for (i=0; i<N; i++) {
        F_re.params = &s[i];
        F_im.params = &s[i];
        pts[1] = s[i];
        
        gsl_integration_qagp(&F_re,pts,3,eps,eps,N_int,w,&temp.re,&error);
        gsl_integration_qagp(&F_im,pts,3,eps,eps,N_int,w,&temp.im,&error);
        
        M.re = gsl_spline_eval(M_etapi_hat.re,s[i],acc);
        M.im = gsl_spline_eval(M_etapi_hat.im,s[i],acc);
        
        log_1 = log(s0/L2*(L2-s[i])/(s[i]-s0));
        log_2 = log(L2/(L2-s[i]));
        
        M_etapi[i].re = temp.re+(-M_PI*M.im+M.re*log_1+Mconst.re*log_2)/s[i];
        M_etapi[i].im = temp.im+(M_PI*M.re+M.im*log_1+Mconst.im*log_2)/s[i];
        
        M_etapi[i].re *= pow(s[i],(double)n)/M_PI;
        M_etapi[i].im *= pow(s[i],(double)n)/M_PI;
    }
    
    gsl_interp_accel_free(acc);
    gsl_integration_workspace_free(w);
}

//etapi_disc_finite function for given s in [s_min,s_etapi] below the cut
void etapi_disc_finite_below_cut(complex *M_etapi, complex_spline M_etapi_hat, double s0, double L2, double s_min, double *s, int N, int n) {
    int i,N_int;
    double error,eps,log_1;
    complex temp,Mconst;
    
    N_int  = 1500;
    eps = 1.0e-14;
    
    gsl_interp_accel *acc = gsl_interp_accel_alloc();
    gsl_integration_workspace *w = gsl_integration_workspace_alloc(N_int);
    gsl_function F_re;
    gsl_function F_im;
    
    double f_re (double z, void * params) {
        double x = *(double *) params;
        double fct = gsl_spline_eval(M_etapi_hat.re,z,acc)/(z*(z-x));
        return fct;
    }
    double f_im (double z, void * params) {
        double x = *(double *) params;
        double fct = gsl_spline_eval(M_etapi_hat.im,z,acc)/(z*(z-x));
        return fct;
    }
    
    F_re.function = &f_re;
    F_im.function = &f_im;
    
    Mconst.re = gsl_spline_eval(M_etapi_hat.re,L2,acc);
    Mconst.im = gsl_spline_eval(M_etapi_hat.im,L2,acc);
    
    for (i=0; i<N; i++) {
        F_re.params = &s[i];
        F_im.params = &s[i];
        
        gsl_integration_qags(&F_re,s0,L2,eps,eps,N_int,w,&temp.re,&error);
        gsl_integration_qags(&F_im,s0,L2,eps,eps,N_int,w,&temp.im,&error);
        
        if (fabs(s[i])<1.e-5) {
            log_1 = 1./L2+s[i]/(2.*L2*L2)+s[i]*s[i]/(3.*L2*L2*L2)+s[i]*s[i]*s[i]/(4.*L2*L2*L2*L2)+s[i]*s[i]*s[i]*s[i]/(5.*L2*L2*L2*L2*L2);
        }
        else {
            log_1 = log(L2/(L2-s[i]))/s[i];
        }
        
        M_etapi[i].re = temp.re+Mconst.re*log_1;
        M_etapi[i].im = temp.im+Mconst.im*log_1;
        
        M_etapi[i].re *= pow(s[i],(double)n)/M_PI;
        M_etapi[i].im *= pow(s[i],(double)n)/M_PI;
    }
    
    gsl_interp_accel_free(acc);
    gsl_integration_workspace_free(w);
}

//etapi_disc_finite function for complex values of s
void etapi_disc_finite_complex(complex *M_etapi, complex_spline M_etapi_hat, double s0, double L2, complex *s, int N, int n) {
    int i,N_int;
    double error,eps,r,phi;
    complex Mconst,temp,log_1;
    
    N_int  = 1500;
    eps = 1.0e-14;
    
    gsl_interp_accel *acc = gsl_interp_accel_alloc();
    gsl_integration_workspace *w = gsl_integration_workspace_alloc(N_int);
    gsl_function F_re,F_im,G_re,G_im;
    
    double f_re (double z, void * params) {
        complex x = *(complex *) params;
        double denom = z*(pow(z-x.re,2.0)+x.im*x.im);
        double fct = gsl_spline_eval(M_etapi_hat.re,z,acc)*(z-x.re)/denom-gsl_spline_eval(M_etapi_hat.im,z,acc)*x.im/denom;
        return fct;
    }
    
    double f_im (double z, void * params) {
        complex x = *(complex *) params;
        double denom = z*(pow(z-x.re,2.0)+x.im*x.im);
        double fct = gsl_spline_eval(M_etapi_hat.re,z,acc)*x.im/denom+gsl_spline_eval(M_etapi_hat.im,z,acc)*(z-x.re)/denom;
        return fct;
    }
    
    Mconst.re = gsl_spline_eval(M_etapi_hat.re,L2,acc);
    Mconst.im = gsl_spline_eval(M_etapi_hat.im,L2,acc);
    
    F_re.function = &f_re;
    F_im.function = &f_im;
    for (i=0; i<N; i++) {
        F_re.params = &s[i];
        F_im.params = &s[i];
        G_re.params = &s[i];
        G_im.params = &s[i];
        
        gsl_integration_qags(&F_re,s0,L2,eps,eps,N_int,w,&temp.re,&error);
        gsl_integration_qags(&F_im,s0,L2,eps,eps,N_int,w,&temp.im,&error);
        
        log_1.re = 0.5*(M_PI*fabs(s[i].im)+2.*s[i].im*atan((s[i].re-L2)/s[i].im)-s[i].re*(log(s[i].im*s[i].im+pow(L2-s[i].re,2.))-2.*log(L2)))/(s[i].re*s[i].re+s[i].im*s[i].im);
        
        log_1.im = 0.5*(M_PI*s[i].re*sign(s[i].im)+2.*s[i].re*atan((s[i].re-L2)/s[i].im)+s[i].im*(log(s[i].im*s[i].im+pow(L2-s[i].re,2.))-2.*log(L2)))/(s[i].re*s[i].re+s[i].im*s[i].im);
        
        M_etapi[i].re = temp.re+log_1.re*Mconst.re-log_1.im*Mconst.im;
        M_etapi[i].im = temp.im+log_1.re*Mconst.im+log_1.im*Mconst.re;
        
        r = pow(sqrt(s[i].re*s[i].re+s[i].im*s[i].im),(double)n)/M_PI;
        phi = n*atan2(s[i].im,s[i].re);
        
        temp.re = r*(M_etapi[i].re*cos(phi)-M_etapi[i].im*sin(phi));
        temp.im = r*(M_etapi[i].re*sin(phi)+M_etapi[i].im*cos(phi));
        
        M_etapi[i].re = temp.re;
        M_etapi[i].im = temp.im;
    }
    
    gsl_interp_accel_free(acc);
    gsl_integration_workspace_free(w);
}

void build_etapi_inhomogeneity(complex **M0_inhom, complex **M1_inhom, complex **M2_inhom, complex_spline M0_finite, complex_spline M2_finite, complex_spline *M0_sing, complex_spline *M2_sing, complex *F0, complex *F2, complex *G0, complex *G2, double *s_cv_plus, double *s_cv_minus, double *s_below_cut, complex *s_complex, double s_min, double L2, int *N, int n0_sub, int n1_sub, int n2_sub) {
    
    int i,j,i_etapi;
    
    complex **M0_etapi = (complex **)malloc(4*sizeof(complex *));
    for (i=0; i<4; i++) {
        M0_etapi[i] = (complex *)malloc(N[i]*sizeof(complex));
    }
    
    complex **M2_etapi = (complex **)malloc(4*sizeof(complex *));
    for (i=0; i<4; i++) {
        M2_etapi[i] = (complex *)malloc(N[i]*sizeof(complex));
    }
    
    complex **M0_etapi_hat = (complex **)malloc(4*sizeof(complex *));
    for (i=0; i<4; i++) {
        M0_etapi_hat[i] = (complex *)malloc(N[i]*sizeof(complex));
    }
    
    complex **M2_etapi_hat = (complex **)malloc(4*sizeof(complex *));
    for (i=0; i<4; i++) {
        M2_etapi_hat[i] = (complex *)malloc(N[i]*sizeof(complex));
    }
    
    complex *R12_cv_plus = (complex *)malloc(N[0]*sizeof(complex));
    complex *Q12_cv_plus = (complex *)malloc(N[0]*sizeof(complex));
    complex *Q12_cv_minus = (complex *)malloc(N[1]*sizeof(complex));
    complex *Q12_complex_f = (complex *)malloc(N[2]*sizeof(complex));
    complex *Q12_complex_g = (complex *)malloc(N[2]*sizeof(complex));
    complex *Q12_below_cut = (complex *)malloc(N[3]*sizeof(complex));
    
    i_etapi = 0;
    for (i=0; s_cv_plus[i]<s_etapi; i++) {
        i_etapi++;
    }
    
    for (i=0; i<N[0]; i++) {
        if (i<i_etapi) {
            Q12_cv_plus[i] = Q_sqrt(s_cv_plus[i],s_etapi,L2);
        }
        else {
            R12_cv_plus[i] = R_sqrt(s_cv_plus[i],s_etapi,L2);
            Q12_cv_plus[i] = Q_sqrt(s_cv_plus[i],s_etapi,L2);
        }
    }
    
    for (i=0; i<N[1]; i++) {
        Q12_cv_minus[i] = Q_sqrt(s_cv_minus[i],s_etapi,L2);
    }
    
    for (i=0; i<N[2]; i++) {
        Q12_complex_f[i] = Q_sqrt_complex_f(s_complex[i],s_etapi);
        Q12_complex_g[i] = Q_sqrt_complex_g(s_complex[i],L2);
    }
    
    for (i=0; i<N[3]; i++) {
        Q12_below_cut[i] = Q_sqrt(s_below_cut[i],s_etapi,L2);
    }
    
    etapi_disc_finite_below_cut(M0_etapi[3],M0_finite,s_etapi,L2,s_min,s_below_cut,N[3],n0_sub);
    etapi_disc_finite_below_cut(M2_etapi[3],M2_finite,s_etapi,L2,s_min,s_below_cut,N[3],n2_sub);
    
    inhomogenity_sqrt_below_cut(M0_etapi_hat[3],M0_sing,F0,G0,s_below_cut,Q12_below_cut,s_etapi,L2,N[3],n0_sub);
    inhomogenity_sqrt_below_cut(M2_etapi_hat[3],M2_sing,F2,G2,s_below_cut,Q12_below_cut,s_etapi,L2,N[3],n2_sub);
    
    etapi_disc_finite_below_cut(M0_etapi[0],M0_finite,s_etapi,L2,s_min,s_cv_plus,i_etapi,n0_sub);
    etapi_disc_finite_below_cut(M2_etapi[0],M2_finite,s_etapi,L2,s_min,s_cv_plus,i_etapi,n2_sub);
    
    inhomogenity_sqrt_below_cut(M0_etapi_hat[0],M0_sing,F0,G0,s_cv_plus,Q12_cv_plus,s_etapi,L2,i_etapi,n0_sub);
    inhomogenity_sqrt_below_cut(M2_etapi_hat[0],M2_sing,F2,G2,s_cv_plus,Q12_cv_plus,s_etapi,L2,i_etapi,n2_sub);
    
    etapi_disc_finite_below_cut(M0_etapi[1],M0_finite,s_etapi,L2,s_min,s_cv_minus,N[1],n0_sub);
    etapi_disc_finite_below_cut(M2_etapi[1],M2_finite,s_etapi,L2,s_min,s_cv_minus,N[1],n2_sub);
    
    inhomogenity_sqrt_below_cut(M0_etapi_hat[1],M0_sing,F0,G0,s_cv_minus,Q12_cv_minus,s_etapi,L2,N[1],n0_sub);
    inhomogenity_sqrt_below_cut(M2_etapi_hat[1],M2_sing,F2,G2,s_cv_minus,Q12_cv_minus,s_etapi,L2,N[1],n2_sub);
    
    etapi_disc_finite_cv_plus(M0_etapi[0]+i_etapi,M0_finite,s_etapi,L2,s_cv_plus+i_etapi,N[0]-i_etapi,n0_sub);
    etapi_disc_finite_cv_plus(M2_etapi[0]+i_etapi,M2_finite,s_etapi,L2,s_cv_plus+i_etapi,N[0]-i_etapi,n2_sub);
    
    inhomogenity_sqrt_cv_plus(M0_etapi_hat[0]+i_etapi,M0_sing,F0,G0,s_cv_plus+i_etapi,R12_cv_plus+i_etapi,Q12_cv_plus+i_etapi,s_etapi,L2,1.,N[0]-i_etapi,n0_sub);
    inhomogenity_sqrt_cv_plus(M2_etapi_hat[0]+i_etapi,M2_sing,F2,G2,s_cv_plus+i_etapi,R12_cv_plus+i_etapi,Q12_cv_plus+i_etapi,s_etapi,L2,1.,N[0]-i_etapi,n2_sub);
    
    etapi_disc_finite_complex(M0_etapi[2],M0_finite,s_etapi,L2,s_complex,N[2],n0_sub);
    etapi_disc_finite_complex(M2_etapi[2],M2_finite,s_etapi,L2,s_complex,N[2],n2_sub);
    
    inhomogenity_sqrt_complex(M0_etapi_hat[2],M0_sing,F0,G0,s_complex,Q12_complex_f,Q12_complex_g,s_etapi,L2,N[2],n0_sub);
    inhomogenity_sqrt_complex(M2_etapi_hat[2],M2_sing,F0,G0,s_complex,Q12_complex_f,Q12_complex_g,s_etapi,L2,N[2],n2_sub);
    
    for (i=0; i<4; i++) {
        for (j=0; j<N[i]; j++) {
            M0_inhom[i][j].re = M0_etapi[i][j].re+M0_etapi_hat[i][j].re;
            M0_inhom[i][j].im = M0_etapi[i][j].im+M0_etapi_hat[i][j].im;
            
            M1_inhom[i][j].re = 0.;
            M1_inhom[i][j].im = 0.;
            
            M2_inhom[i][j].re = M2_etapi[i][j].re+M2_etapi_hat[i][j].re;
            M2_inhom[i][j].im = M2_etapi[i][j].im+M2_etapi_hat[i][j].im;
            
            if (isnan(M0_inhom[i][j].re)==1 || isnan(M0_inhom[i][j].im)==1) {
                printf("FATAL ERROR: M0_inhom[%d].re=%g, M0_imhom[%d].im=%g, step=%d\n",i,M0_inhom[i][j].re,i,M0_inhom[i][j].im,j);
                exit(1);
            }
            
            if (isnan(M1_inhom[i][j].re)==1 || isnan(M1_inhom[i][j].im)==1) {
                printf("FATAL ERROR: M1_inhom[%d].re=%g, M1_imhom[%d].im=%g, step=%d\n",i,M1_inhom[i][j].re,i,M1_inhom[i][j].im,j);
                exit(1);
            }
            
            if (isnan(M2_inhom[i][j].re)==1 || isnan(M2_inhom[i][j].im)==1) {
                printf("FATAL ERROR: M2_inhom[%d].re=%g, M2_imhom[%d].im=%g, step=%d\n",i,M2_inhom[i][j].re,i,M2_inhom[i][j].im,j);
                exit(1);
            }
        }
    }
    
    for (i=0; i<4; i++) {
        free(M0_etapi[i]);
        free(M2_etapi[i]);
        free(M0_etapi_hat[i]);
        free(M2_etapi_hat[i]);
    }
    free(M0_etapi);
    free(M2_etapi);
    free(M0_etapi_hat);
    free(M2_etapi_hat);
    
    free(R12_cv_plus);
    free(Q12_cv_plus);
    free(Q12_cv_minus);
    free(Q12_complex_f);
    free(Q12_complex_g);
    free(Q12_below_cut);
}


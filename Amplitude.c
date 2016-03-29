#include "Basic.h"
#include "Amplitude.h"

void subtraction_constant(char *sub_const, int *n0, int *n1, int *n2) {
    if (str_eq_str(sub_const,"a0")==1) {
        *n0 = 0;
        *n1 = -1;
        *n2 = -1;
    }
    else if (str_eq_str(sub_const,"b0")==1) {
        *n0 = 1;
        *n1 = -1;
        *n2 = -1;
    }
    else if (str_eq_str(sub_const,"c0")==1) {
        *n0 = 2;
        *n1 = -1;
        *n2 = -1;
    }
    else if (str_eq_str(sub_const,"d0")==1) {
        *n0 = 3;
        *n1 = -1;
        *n2 = -1;
    }
    else if (str_eq_str(sub_const,"a1")==1) {
        *n0 = -1;
        *n1 = 0;
        *n2 = -1;
    }
    else if (str_eq_str(sub_const,"a2")==1) {
        *n0 = -1;
        *n1 = -1;
        *n2 = 0;
    }
    else if (str_eq_str(sub_const,"b2")==1) {
        *n0 = -1;
        *n1 = -1;
        *n2 = 1;
    }
    else if (str_eq_str(sub_const,"etapi")==1) {
        *n0 = -1;
        *n1 = -1;
        *n2 = -1;
    }
    else {
        printf("FATAL ERROR: Bad input, unknown subtraction constant: %s\n",sub_const);
        exit(1);
    }
}

void build_inhomogenity_init(complex **M_inhom, int *N) {
    int i,j;
    
    for (i=0; i<4; i++) {
        for (j=0; j<N[i]; j++) {
            M_inhom[i][j].re = 0.0;
            M_inhom[i][j].im = 0.0;
        }
    }
}

void build_amplitude(complex_spline *M, complex **omnes, complex **M_inhom, double *s_cv_plus, double *s_cv_minus, double *s_below_cut, complex *s_complex, double *phi, double s0, double L2, double s_min, int *N, int n) {
    int i;
    double si,sf;
    complex s,m,O;
    double *M_re = (double *)malloc((N[0]+N[1]+N[2]+N[3])*sizeof(double));
    double *M_im = (double *)malloc((N[0]+N[1]+N[2]+N[3])*sizeof(double));
    
    for (i=0; i<N[0]; i++) {
        O.re = omnes[0][i].re;
        O.im = omnes[0][i].im;
        m.re = M_inhom[0][i].re;
        m.im = M_inhom[0][i].im;
        
        if (isnan(m.re)==1) {
            printf("FATAL ERROR: In function build_amplitude, region 0: m.re=nan, step=%d\n",i);
            exit(1);
        }
        else if (isnan(m.im)==1) {
            printf("FATAL ERROR: In function build_amplitude, region 0: m.im=nan, step=%d\n",i);
            exit(1);
        }
        
        if (isnan(O.re)==1) {
            printf("FATAL ERROR: In function build_amplitude, region 0: O.re=nan, step=%d\n",i);
            exit(1);
        }
        else if (isnan(O.im)==1) {
            printf("FATAL ERROR: In function build_amplitude, region 0: O.im=nan, step=%d\n",i);
            exit(1);
        }
        
        if (n==-1) {
            M_re[i] = O.re*m.re-O.im*m.im;
            M_im[i] = O.re*m.im+O.im*m.re;
        }
        
        else if (n==0) {
            M_re[i] = O.re*(1.0+m.re)-O.im*m.im;
            M_im[i] = O.re*m.im+O.im*(1.0+m.re);
        }
        
        else if (n==1) {
            s.re = s_cv_plus[i];
            
            M_re[i] = O.re*(s.re+m.re)-O.im*m.im;
            M_im[i] = O.re*m.im+O.im*(s.re+m.re);
        }
        
        else if (n==2) {
            s.re = s_cv_plus[i]*s_cv_plus[i];
            
            M_re[i] = O.re*(s.re+m.re)-O.im*m.im;
            M_im[i] = O.re*m.im+O.im*(s.re+m.re);
        }
        
        else if (n==3) {
            s.re = s_cv_plus[i]*s_cv_plus[i]*s_cv_plus[i];
            
            M_re[i] = O.re*(s.re+m.re)-O.im*m.im;
            M_im[i] = O.re*m.im+O.im*(s.re+m.re);
        }
        
    }
    
    si = s_cv_plus[0];
    sf = s_cv_plus[N[0]];
    s_cv_plus[0] = s0;
    s_cv_plus[N[0]] = L2;
    gsl_spline_init(M[0].re,s_cv_plus,M_re,N[0]);
    gsl_spline_init(M[0].im,s_cv_plus,M_im,N[0]);
    s_cv_plus[0] = si;
    s_cv_plus[N[0]] = sf;
    
    for (i=0; i<N[1]; i++) {
        O.re = omnes[1][i].re;
        O.im = omnes[1][i].im;
        m.re = M_inhom[1][i].re;
        m.im = M_inhom[1][i].im;
        
        if (isnan(m.re)==1) {
            printf("FATAL ERROR: In function build_amplitude, region 1: m.re=nan, step=%d\n",i);
            int j;
            for (j=0; j<N[1]; j++) {
                printf("%.3e %.3e %.3e\n",s_cv_minus[j],M_inhom[1][j].re,M_inhom[1][j].im);
            }
            exit(1);
        }
        else if (isnan(m.im)==1) {
            printf("FATAL ERROR: In function build_amplitude, region 1: m.im=nan, step=%d\n",i);
            int j;
            for (j=0; j<N[1]; j++) {
                printf("%.3e %.3e %.3e\n",s_cv_minus[j],M_inhom[1][j].re,M_inhom[1][j].im);
            }
            exit(1);
        }
        
        if (isnan(O.re)==1) {
            printf("FATAL ERROR: In function build_amplitude, region 1: O.re=nan, step=%d\n",i);
            int j;
            for (j=0; j<N[1]; j++) {
                printf("%.3e %.3e %.3e\n",s_cv_minus[j],M_inhom[1][j].re,M_inhom[1][j].im);
            }
            exit(1);
        }
        else if (isnan(O.im)==1) {
            printf("FATAL ERROR: In function build_amplitude, region 1: O.im=nan, step=%d\n",i);
            int j;
            for (j=0; j<N[1]; j++) {
                printf("%.3e %.3e %.3e\n",s_cv_minus[j],M_inhom[1][j].re,M_inhom[1][j].im);
            }
            exit(1);
        }
        
        if (n==-1) {
            M_re[i] = O.re*m.re-O.im*m.im;
            M_im[i] = O.re*m.im+O.im*m.re;
        }
        
        else if (n==0) {
            M_re[i] = O.re*(1.0+m.re)-O.im*m.im;
            M_im[i] = O.re*m.im+O.im*(1.0+m.re);
        }
        
        else if (n==1) {
            s.re = s_cv_minus[i];
            
            M_re[i] = O.re*(s.re+m.re)-O.im*m.im;
            M_im[i] = O.re*m.im+O.im*(s.re+m.re);
        }
        
        else if (n==2) {
            s.re = s_cv_minus[i]*s_cv_minus[i];
            
            M_re[i] = O.re*(s.re+m.re)-O.im*m.im;
            M_im[i] = O.re*m.im+O.im*(s.re+m.re);
        }
        
        else if (n==3) {
            s.re = s_cv_minus[i]*s_cv_minus[i]*s_cv_minus[i];
            
            M_re[i] = O.re*(s.re+m.re)-O.im*m.im;
            M_im[i] = O.re*m.im+O.im*(s.re+m.re);
        }
        
    }
    
    si = s_cv_minus[0];
    sf = s_cv_minus[N[1]];
    s_cv_minus[0] = s0;
    s_cv_minus[N[1]] = (METAP/MPION+1.0);
    gsl_spline_init(M[1].re,s_cv_minus,M_re,N[1]);
    gsl_spline_init(M[1].im,s_cv_minus,M_im,N[1]);
    s_cv_minus[0] = si;
    s_cv_minus[N[1]] = sf;
    
    for (i=0; i<N[2]; i++) {
        O.re = omnes[2][i].re;
        O.im = omnes[2][i].im;
        m.re = M_inhom[2][i].re;
        m.im = M_inhom[2][i].im;
        
        if (isnan(m.re)==1) {
            printf("FATAL ERROR: In function build_amplitude, region 2: m.re=nan, step=%d\n",i);
            exit(1);
        }
        else if (isnan(m.im)==1) {
            printf("FATAL ERROR: In function build_amplitude, region 2: m.im=nan, step=%d\n",i);
            exit(1);
        }
        
        if (isnan(O.re)==1) {
            printf("FATAL ERROR: In function build_amplitude, region 2: O.re=nan, step=%d\n",i);
            exit(1);
        }
        else if (isnan(O.im)==1) {
            printf("FATAL ERROR: In function build_amplitude, region 2: O.im=nan, step=%d\n",i);
            exit(1);
        }
        
        if (n==-1) {
            M_re[i] = O.re*m.re-O.im*m.im;
            M_im[i] = O.re*m.im+O.im*m.re;
        }
        
        else if (n==0) {
            M_re[i] = O.re*(1.0+m.re)-O.im*m.im;
            M_im[i] = O.re*m.im+O.im*(1.0+m.re);
        }
        
        else if (n==1) {
            s = s_complex[i];
            
            M_re[i] = O.re*(s.re+m.re)-O.im*(s.im+m.im);
            M_im[i] = O.re*(s.im+m.im)+O.im*(s.re+m.re);
        }
        
        else if (n==2) {
            s.re = s_complex[i].re*s_complex[i].re-s_complex[i].im*s_complex[i].im;
            s.im = 2.0*s_complex[i].re*s_complex[i].im;
            
            M_re[i] = O.re*(s.re+m.re)-O.im*(s.im+m.im);
            M_im[i] = O.re*(s.im+m.im)+O.im*(s.re+m.re);
        }
        
        else if (n==3) {
            s.re = s_complex[i].re*s_complex[i].re*s_complex[i].re-3.0*s_complex[i].re*s_complex[i].im*s_complex[i].im;
            s.im = 3.0*s_complex[i].re*s_complex[i].re*s_complex[i].im-s_complex[i].im*s_complex[i].im*s_complex[i].im;
            
            M_re[i] = O.re*(s.re+m.re)-O.im*(s.im+m.im);
            M_im[i] = O.re*(s.im+m.im)+O.im*(s.re+m.re);
        }
    }
    
    si = phi[0];
    sf = phi[N[2]];
    phi[0] = 0.0;
    phi[N[2]] = 2.0*M_PI;
    gsl_spline_init(M[2].re,phi,M_re,N[2]);
    gsl_spline_init(M[2].im,phi,M_im,N[2]);
    phi[0] = si;
    phi[N[2]] = sf;
    
    for (i=0; i<N[3]; i++) {
        O.re = omnes[3][i].re;
        O.im = omnes[3][i].im;
        m.re = M_inhom[3][i].re;
        m.im = M_inhom[3][i].im;
        
        if (isnan(m.re)==1) {
            printf("FATAL ERROR: In function build_amplitude, region 3: m.re=nan, step=%d\n",i);
            exit(1);
        }
        else if (isnan(m.im)==1) {
            printf("FATAL ERROR: In function build_amplitude, region 3: m.im=nan, step=%d\n",i);
            exit(1);
        }
        
        if (isnan(O.re)==1) {
            printf("FATAL ERROR: In function build_amplitude, region 3: O.re=nan, step=%d\n",i);
            exit(1);
        }
        else if (isnan(O.im)==1) {
            printf("FATAL ERROR: In function build_amplitude, region 3: O.im=nan, step=%d\n",i);
            exit(1);
        }
        
        if (n==-1) {
            M_re[i] = O.re*m.re-O.im*m.im;
            M_im[i] = O.re*m.im+O.im*m.re;
        }
        
        else if (n==0) {
            M_re[i] = O.re*(1.0+m.re)-O.im*m.im;
            M_im[i] = O.re*m.im+O.im*(1.0+m.re);
        }
        
        else if (n==1) {
            s.re = s_below_cut[i];
            
            M_re[i] = O.re*(s.re+m.re)-O.im*m.im;
            M_im[i] = O.re*m.im+O.im*(s.re+m.re);
        }
        
        else if (n==2) {
            s.re = s_below_cut[i]*s_below_cut[i];
            
            M_re[i] = O.re*(s.re+m.re)-O.im*m.im;
            M_im[i] = O.re*m.im+O.im*(s.re+m.re);
        }
        
        else if (n==3) {
            s.re = s_below_cut[i]*s_below_cut[i]*s_below_cut[i];
            
            M_re[i] = O.re*(s.re+m.re)-O.im*m.im;
            M_im[i] = O.re*m.im+O.im*(s.re+m.re);
        }
    }
    
    si = s_below_cut[0];
    sf = s_below_cut[N[3]];
    s_below_cut[0] = s_min;
    s_below_cut[N[3]] = s0;
    gsl_spline_init(M[3].re,s_below_cut,M_re,N[3]);
    gsl_spline_init(M[3].im,s_below_cut,M_im,N[3]);
    s_below_cut[0] = si;
    s_below_cut[N[3]] = sf;
    
    free(M_re);
    free(M_im);
}

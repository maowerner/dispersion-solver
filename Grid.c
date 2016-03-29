#include "Basic.h"
#include "Grid.h"

void grid_point_size(int *N, int Nsin, double s_step, double s_step_low_high, double s_min, double s_max, double s_th, double eps_b) {
    
    //below cut
    //points from s_min to 0
    N[0] = (int)ceil(-s_min/s_step_low_high);
    
    //points from 0 to s_th
    N[1] = (int)ceil(s_th/s_step)+Nsin;
    
    //above cut cv plus
    //points from s_th to 0.5*(METAP^2-MPION^2)
    N[2] = (int)ceil((0.5*(METAP_SQUARED-1.0)-s_th)/s_step)+Nsin;
    
    if (METAP!=META) {
        //points from 0.5*(METAP^2-MPION^2) to (META+MPION)^2
        N[3] = (int)ceil((s_etapi-0.5*(METAP_SQUARED-1.0))/s_step)+Nsin;
        
        //points from (META+MPION)^2 to (METAP-MPION)^2
        N[9] = (int)ceil((METAP_M_MPION_SQUARED-s_etapi)/s_step)+2*Nsin;
    }
    else {
        //points from 0.5*(METAP^2-MPION^2) to (METAP-MPION)^2
        N[3] = (int)ceil((METAP_M_MPION_SQUARED-0.5*(METAP_SQUARED-1.0))/s_step)+Nsin;
        N[9] = 0;
    }
    
    //points from (METAP-MPION)^2 to (METAP+MPION)^2
    N[4] = (int)ceil(((METAP_P_MPION_SQUARED-eps_b)-METAP_M_MPION_SQUARED)/s_step)+Nsin;
    
    //points from (METAP+MPION)^2 to 100*MPION^2
    N[5] = (int)ceil((100.-(METAP_P_MPION_SQUARED+eps_b))/s_step);
    
    //points from 100*MPION^2 to s_max
    N[6] = (int)ceil((s_max-100.)/s_step_low_high);
    
    //above cut cv minus
    //points from s_th to MPION*(METAP+MPION)
    N[7] = (int)ceil((METAP/MPION+1.0-s_th)/s_step)+2*Nsin;
    
    //complex plane
    N[8] = (int)ceil(integration_III_circumference()/s_step)+2*Nsin;
}

void build_s_below_cut(double *s, double s_min, double s_th, double eps_a, int *N, int Nsin) {
    int i,n,m;
    double si,sf,ds;
    
    n = N[0];
    si = s_min;
    sf = 0.;
    ds = (sf-si)/n;
    s[0] = si+eps_a;
    for (i=1; i<n; i++) {
        s[i] = si+ds*i;
    }
    s[i] = sf;
    m = N[1]-Nsin;
    si = 0.;
    sf = s_th;
    ds = (sf-si)/m;
    for (i=1; i<m; i++) {
        s[i+n] = si+ds*i;
    }
    for (i=0; i<Nsin; i++) {
        s[n+m+i] = s[n+m+i-1]+ds*pow(0.5,i+1.);
    }
}

void build_s_cv_plus(double *s, double s_th, double s_max, double eps_a, double eps_b, int *N, int Nsin) {
    int i,n,m;
    double si,sf,ds;
    
    n = N[2]-Nsin;
    si = s_th;
    sf = 0.5*(METAP_SQUARED-1.0);
    ds = (sf-si)/n;
    s[0] = si+ds*pow(0.5,Nsin+1.);
    for (i=1; i<Nsin; i++) {
        s[i] = s[i-1]+ds*pow(0.5,Nsin+1.-i);
    }
    for (i=0; i<n-1; i++) {
        s[i+Nsin] = si+ds*(i+1.);
    }
    s[N[2]-1] = sf-eps_a;
    if (METAP!=META) {
        n = N[3]-Nsin;
        m = N[2];
        si = 0.5*(METAP_SQUARED-1.0);
        sf = s_etapi;
        ds = (sf-si)/n;
        s[m] = si+eps_a;
        for (i=1; i<n; i++) {
            s[i+m] = si+ds*i;
        }
        m += n;
        for (i=0; i<Nsin; i++) {
            s[m+i] = s[m+i-1]+ds*pow(0.5,i+1.);
        }
        n = N[9]-2*Nsin;
        m = N[2]+N[3];
        si = s_etapi;
        sf = METAP_M_MPION_SQUARED;
        ds = (sf-si)/(n+1.);
        s[m] = si+ds*pow(0.5,Nsin+1.);
        for (i=1; i<Nsin; i++) {
            s[m+i] = s[m+i-1]+ds*pow(0.5,Nsin+1.-i);
        }
        m += Nsin;
        for (i=0; i<n; i++) {
            s[m+i] = si+ds*(i+1.);
        }
        m += n;
        for (i=0; i<Nsin; i++) {
            s[m+i] = s[m+i-1]+ds*pow(0.5,i+1.);
        }
    }
    else {
        n = N[3]-Nsin;
        m = N[2];
        si = 0.5*(METAP_SQUARED-1.0);
        sf = METAP_M_MPION_SQUARED;
        ds = (sf-si)/n;
        s[m] = si+eps_a;
        for (i=1; i<n; i++) {
            s[i+m] = si+ds*i;
        }
        for (i=0; i<Nsin; i++) {
            s[n+m+i] = s[n+m+i-1]+ds*pow(0.5,i+1.);
        }
    }
    n = N[4]-Nsin;
    m = N[2]+N[3]+N[9];
    si = METAP_M_MPION_SQUARED;
    sf = METAP_P_MPION_SQUARED-eps_b;
    ds = (sf-si)/n;
    s[m] = si+ds*pow(0.5,Nsin+1.);
    for (i=1; i<Nsin; i++) {
        s[m+i] = s[m+i-1]+ds*pow(0.5,Nsin+1.-i);
    }
    m += Nsin;
    for (i=0; i<n; i++) {
        s[m+i] = si+ds*(i+1.);
    }
    n = N[5];
    m = N[2]+N[3]+N[9]+N[4];
    si = METAP_P_MPION_SQUARED+eps_b;
    sf = 100.;
    ds = (sf-si)/(n-1.);
    for (i=0; i<n; i++) {
        s[m+i] = si+ds*i;
    }
    n = N[6];
    m = N[2]+N[3]+N[9]+N[4]+N[5];
    si = 100.;
    sf = s_max;
    ds = (sf-si)/n;
    for (i=0; i<n-1; i++) {
        s[m+i] = si+ds*(i+1.);
    }
    s[m+n-1] = s_max-eps_a;
}

void build_s_cv_minus(double *s, double s_th, int *N, int Nsin) {
    int i,n,m;
    double si,sf,ds;
    
    n = N[7]-2*Nsin;
    si = s_th;
    sf = METAP/MPION+1.0;
    ds = (sf-si)/(n+1.);
    s[0] = si+ds*pow(0.5,Nsin+1.);
    for (i=1; i<Nsin; i++) {
        s[i] = s[i-1]+ds*pow(0.5,Nsin+1.-i);
    }
    m = Nsin;
    for (i=0; i<n; i++) {
        s[m+i] = si+ds*(i+1.);
    }
    m += n;
    for (i=0; i<Nsin; i++) {
        s[m+i] = s[m+i-1]+ds*pow(0.5,i+1.);
    }
    s[N[7]-1] = sf;
}

void build_s_complex(complex *s, double *phi, int *N, int Nsin) {
    int i,n,m;
    double dphi;
    
    n = N[8]-2*Nsin;
    dphi = 2.0*M_PI/(n+1.);
    phi[0] = dphi*pow(0.5,Nsin+1.);
    s[0] = integration_III_gamma_phi(phi[0]);
    for (i=1; i<Nsin; i++) {
        phi[i] = phi[i-1]+dphi*pow(0.5,Nsin+1.-i);
        s[i] = integration_III_gamma_phi(phi[i]);
    }
    m = Nsin;
    for (i=0; i<n; i++) {
        phi[m+i] = dphi*(i+1.);
        s[m+i] = integration_III_gamma_phi(phi[m+i]);
    }
    m += n;
    for (i=0; i<Nsin; i++) {
        phi[m+i] = phi[m+i-1]+dphi*pow(0.5,i+1.);
        s[m+i] = integration_III_gamma_phi(phi[m+i]);
    }
    
    phi[0] = 0.;
    phi[N[8]-1] = 2.*M_PI;
}

void build_grid(double *s_cv_plus, double *s_cv_minus, double *s_below_cut, complex *s_complex, double *phi, double s_min, double s_max, double s_th, double eps_a, double eps_b, int *N, int Nsin) {
    int i,n;
    
    build_s_below_cut(s_below_cut,s_min,s_th,eps_a,N,Nsin);
    n = N[0]+N[1];
    for (i=0; i<n; i++) {
        if (s_below_cut[i]>=s_below_cut[i+1] && i+1<n-1) {
            printf("FATAL ERROR: s_below_cut is not stictly monotonically increasing! step=%d, s=%.10e\n\n",i,s_below_cut[i]);
            exit(1);
        }
    }
    
    build_s_cv_plus(s_cv_plus,s_th,s_max,eps_a,eps_b,N,Nsin);
    n = N[2]+N[3]+N[4]+N[5]+N[6]+N[9];
    for (i=0; i<n; i++) {
        if (s_cv_plus[i]>=s_cv_plus[i+1] && i+1<n-1) {
            printf("FATAL ERROR: s_cv_plus is not stictly monotonically increasing! step=%d, s=%.10e\n\n",i,s_cv_plus[i]);
            exit(1);
        }
    }
    
    build_s_cv_minus(s_cv_minus,s_th,N,Nsin);
    n = N[7];
    for (i=0; i<n; i++) {
        if (s_cv_minus[i]>=s_cv_minus[i+1] && i+1<n-1) {
            printf("FATAL ERROR: s_cv_minus is not stictly monotonically increasing! step=%d, s=%.10e\n\n",i,s_cv_minus[i]);
            exit(1);
        }
    }
    
    build_s_complex(s_complex,phi,N,Nsin);
    n = N[8];
    for (i=0; i<n; i++) {
        if (phi[i]>=phi[i+1] && i+1<n-1) {
            printf("FATAL ERROR: phi is not stictly monotonically increasing! step=%d, phi=%.10e\n\n",i,phi[i]);
            exit(1);
        }
    }
}

//
void build_s_etapi_disc(double *s_etapi_disc, double s_init, double s_final, double eps, int *N, int *N_temp, int Nsin) {
    int i,n,m;
    double s,ds;
    
    s_etapi_disc[0] = s_init+eps;
    n = N[0]-Nsin;
    ds = (METAP_M_MPION_SQUARED-s_init)/n;
    for (i=1; i<n; i++) {
        s_etapi_disc[i] = s_init+ds*i;
    }
    for (i=0; i<Nsin; i++) {
        s_etapi_disc[n+i] = s_etapi_disc[n-1+i]+ds*pow(0.5,i+1);
    }
    
    m = N[0]+Nsin;
    n = N_temp[0]-Nsin;
    ds = (100.-METAP_M_MPION_SQUARED)/n;
    for (i=0; i<n; i++) {
        s_etapi_disc[m+i] = METAP_M_MPION_SQUARED+ds*(i+1);
    }
    m = N[0];
    for (i=Nsin-1; i>=0; i--) {
        s_etapi_disc[m+i] = s_etapi_disc[m+i+1]-ds*pow(0.5,Nsin-i);
    }
    n = N_temp[1];
    ds = (s_final-100.)/n;
    m = N[0]+N_temp[0];
    for (i=0; i<n; i++) {
        s_etapi_disc[m+i] = 100.+ds*(i+1);
    }
    
    s_etapi_disc[N[0]+N[1]-1] = s_final-eps;
    //    printf("%.10e\n%.10e\n%.10e\n",s_etapi_disc[N[0]-1],METAP_M_MPION_SQUARED,s_etapi_disc[N[0]]);
    //    for (i=0; i<(N[0]+N[1]); i++) {
    //        printf("%d %+.10e\n",i,s_etapi_disc[i]);
    //        if (s_etapi_disc[i]>s_etapi_disc[i+1]) {
    //            printf("%d nan %d\n",i-N[0],i-N[0]-Nsin);
    //        }
    //    }
}

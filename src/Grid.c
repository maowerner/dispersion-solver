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
    N[8] = 100;//(int)ceil(integration_III_circumference()/s_step)+2*Nsin;
}

void Ngrid_etap_eta2pi(int *N, int Nsin, double s_step, double s_step_low_high, double s_min, double s_max, double s_th, double a_s, double b_s, double s_ps, double s_m, double s_12, double t_th, double a_t, double b_t, double t_ps, double t_m, double t_12, double u_m, double u_12, double eps_b) {
    
    //left from cut
    //points from s_min to 0
    N[0] = (int)ceil(-s_min/s_step_low_high);
    
    //points from 0 to s_th(pi0pi0)
    N[1] = (int)ceil(4.*pow(MPION0/MPION,2.)/s_step)+Nsin;
    
    //points from s_th(pi0pi0) to s_th(pi+pi-)
    N[2] = (int)ceil((4.-4.*pow(MPION0/MPION,2.))/(0.1*s_step))+2*Nsin;
    
    //on the upper rim of the cut
    //points from s_th(pi+pi-) to s_12
    N[3] = (int)ceil((s_12-4.)/s_step)+Nsin;
    
    //points from s_12 to a_s
    N[4] = (int)ceil((a_s-s_12)/s_step)+Nsin;
    
    //points from a_s to t_th
    N[5] = (int)ceil((t_th-a_s)/s_step)+2*Nsin;
    
    //points from t_th to t_12
    N[6] = (int)ceil((t_12-t_th)/s_step)+Nsin;
    
    //points from t_12 to u_12
    N[7] = (int)ceil((u_12-t_12)/s_step);
    
    //points from u_12 to a_t
    N[8] = (int)ceil((a_t-u_12)/s_step)+Nsin;
    
    //points from a_t to KKbar cusp
    N[9] = (int)ceil((49.-a_t)/s_step)+Nsin;
    
    //points around the KKbar cusp [49,51]MPION^2
    N[10] = (int)ceil((51.-49.)/(0.01*s_step));
    
    //points from KKbar cusp to b_t-eps
    N[11] = (int)ceil(((b_t-eps_b)-51.)/s_step)+Nsin;
    
    //points from b_t+eps_b to b_s-eps_b
    N[12] = (int)ceil(((b_s-eps_b)-(b_t+eps_b))/s_step);
    
    //points from b_s+eps_b to s_max
    N[13] = (int)ceil((s_max-(b_s+eps_b))/s_step_low_high);
    
    //on the lower rim of the cut
    //points from s_th to t_m
    N[14] = (int)ceil((t_m-s_th)/s_step)+Nsin;
    
    //points from t_th to s_m
    N[15] = (int)ceil((s_m-t_th)/s_step)+Nsin;
    
    //complex plane for s-channel sp-variable
    N[16] = (int)ceil(M_PI*fabs(cmp_im(gamma_path(0.5*(a_t+b_t),a_t,b_t,t_ps,t_th,'-','-')))/(2.*s_step));
    
    //complex plane for t-channel tp-variable
    N[17] = (int)ceil(M_PI*fabs(cmp_im(gamma_path(0.5*(a_s+b_s),a_s,b_s,s_ps,s_th,'0','-')))/(2.*s_step));
    
    //complex plane for t-channel up-variable
    N[18] = (int)ceil(M_PI*fabs(cmp_im(gamma_path(0.5*(a_t+b_t),a_t,b_t,t_ps,t_th,'+','-')))/(2.*s_step));;
}

//build interval of the form [si,sf]
void closed_to_closed(double *s, double si, double sf, int N) {
    int i;
    double ds = (sf-si)/(N-1.);
    
    s[0] = si;
    for (i=1; i<N-1; i++) {
        s[i] = si+ds*i;
    }
    s[i] = sf;
}

//build interval of the form [si,sf[ where Nsin points approach sf by (1/2)^Nsin the distance to sf
void closed_to_open(double *s, double si, double sf, int N, int Nsin) {
    int i,n;
    double ds;
    
    n = N-Nsin;
    ds = (sf-si)/n;
    
    s[0] = si;
    for (i=1; i<n; i++) {
        s[i] = si+ds*i;
    }
    for (i=0; i<Nsin; i++) {
        s[n+i] = s[n+i-1]+ds*pow(0.5,i+1.);
    }
}

//build interval of the form ]si,sf] where Nsin points approach si by (1/2)^Nsin the distance to si
void open_to_closed(double *s, double si, double sf, int N, int Nsin) {
    int i,n;
    double ds;
    
    n = N-Nsin;
    ds = (sf-si)/n;
    s[0] = si+ds*pow(0.5,Nsin+1.);
    for (i=1; i<Nsin; i++) {
        s[i] = s[i-1]+ds*pow(0.5,Nsin+1.-i);
    }
    for (i=0; i<n-1; i++) {
        s[i+Nsin] = si+ds*(i+1.);
    }
    s[N-1] = sf;
}

//build interval of the form ]si,sf[ where Nsin points approach si/sf by (1/2)^Nsin the distance to si/sf
void open_to_open(double *s, double si, double sf, int N, int Nsin) {
    int i,n;
    double ds;
    
    n = N-2*Nsin;
    ds = (sf-si)/(n+1.);
    s[0] = si+ds*pow(0.5,Nsin+1.);
    for (i=1; i<Nsin; i++) {
        s[i] = s[i-1]+ds*pow(0.5,Nsin+1.-i);
    }
    for (i=0; i<n; i++) {
        s[Nsin+i] = si+ds*(i+1.);
    }
    for (i=0; i<Nsin; i++) {
        s[Nsin+n+i] = s[Nsin+n+i-1]+ds*pow(0.5,i+1.);
    }
}

void grid_etap_eta2pi(double **s, double **ys, double **yt, double **yu, complex **s_cmp, complex **t_cmp, complex **u_cmp, int *N, int Nsin, double s_step, double s_step_low_high, double s_min, double s_max, double s_th, double a_s, double b_s, double s_ps, double s_m, double s_12, double t_th, double a_t, double b_t, double t_ps, double t_m, double t_12, double u_m, double u_12, double eps_b, double eps_a) {
    int i;
    double si,sf,ds,y_min,y_max,eps;
    
    eps = 1.e-10;
    
    //s in range [smin,0]
    closed_to_closed(s[0],s_min,-eps_a,N[0]);
    
    //s in range [0,s_th(pi0pi0)]
    closed_to_open(s[1],eps_a,4.*pow(MPION0/MPION,2.),N[1],Nsin);
    
    //s in range [s_th(pi0pi0),s_th(pi+pi-)]
    open_to_open(s[2],4.*pow(MPION0/MPION,2.),4.,N[2],Nsin);
    
    //s in range [s_th,s_12]
    open_to_closed(s[3],4.,s_12-eps_a,N[3],Nsin);
    
    //s in range [s_12,a_s]
    closed_to_open(s[4],s_12+eps_a,a_s,N[4],Nsin);
    
    //s in range [a_s,t_th]
    open_to_open(s[5],a_s,t_th,N[5],Nsin);
    
    //s in range [t_th,t_12]
    open_to_closed(s[6],t_th,t_12-eps_a,N[6],Nsin);
    
    //s in range [t_12,u_12]
    closed_to_closed(s[7],t_12+eps_a,u_12-eps_a,N[7]);
    
    //s in range [u_12,a_t]
    closed_to_open(s[8],u_12+eps_a,a_t,N[8],Nsin);
    
    //s in range [a_t,KKbar cusp]
    open_to_closed(s[9],a_t,49.-eps_a,N[9],Nsin);
    
    //s in range [KKbar cusp-1,KKbar cusp+1]
    closed_to_closed(s[10],49.+eps_a,51.-eps_a,N[10]);
    
    //s in range [KKbar cusp,b_t]
    closed_to_closed(s[11],51.+eps_a,b_t-eps_b,N[11]);
    
    //s in range [b_t,b_s]
    closed_to_closed(s[12],b_t+eps_b,b_s-eps_b,N[12]);
    
    //s in range [b_s,s_max]
    closed_to_closed(s[13],b_s+eps_b,s_max-eps_a,N[13]);
    
    //s in range [s_th,t_m]
    open_to_closed(s[14],s_th,t_m,N[14],Nsin);
    
    //s in range [t_th,s_m]
    open_to_closed(s[15],t_th,s_m,N[15],Nsin);
    
    //complex plane for s-channel sp-variable
    //ya
    y_min = 0.5*(sqrt(sqrt(2.*(a_t+b_t))+2.*sqrt(a_t))-sqrt(sqrt(2.*(a_t+b_t))-2.*sqrt(a_t)));
    y_max = pow(a_t,1./4.);
    closed_to_closed(ys[0],y_min,y_max,N[16]);
    //yb
    y_min = 0.5*(sqrt(2.*sqrt(b_t)+sqrt(2.*(a_t+b_t)))-sqrt(2.*sqrt(b_t)-sqrt(2.*(a_t+b_t))));
    y_max = pow(b_t,1./4.);
    closed_to_closed(ys[1],y_min,y_max,N[16]);
    
    for (i=0; i<N[16]; i++) {
        s_cmp[0][i] = gamma_III(ys[0][i],a_t,b_t,t_ps,t_th,'a','-','-');
        s_cmp[1][i] = gamma_III(ys[1][i],a_t,b_t,t_ps,t_th,'b','-','-');
    }
    s_cmp[0][N[16]-1].im = -eps;
    s_cmp[1][N[16]-1].im = -eps;
    
    //complex plane for t-channel tp-variable
    //ya
    y_min = 0.5*(sqrt(sqrt(2.*(a_s+b_s))+2.*sqrt(a_s))-sqrt(sqrt(2.*(a_s+b_s))-2.*sqrt(a_s)));
    y_max = pow(a_s,1./4.);
    closed_to_closed(yt[0],y_min,y_max,N[17]);
    //yb
    y_min = 0.5*(sqrt(2.*sqrt(b_s)+sqrt(2.*(a_s+b_s)))-sqrt(2.*sqrt(b_s)-sqrt(2.*(a_s+b_s))));
    y_max = pow(b_s,1./4.);
    closed_to_closed(yt[1],y_min,y_max,N[17]);
    
    for (i=0; i<N[17]; i++) {
        t_cmp[0][i] = gamma_III(yt[0][i],a_s,b_s,s_ps,s_th,'a','0','+');
        t_cmp[1][i] = gamma_III(yt[1][i],a_s,b_s,s_ps,s_th,'b','0','+');
    }
    t_cmp[0][N[17]-1].im = eps;
    t_cmp[1][N[17]-1].im = eps;
    
    //complex plane for t-channel up-variable
    //ya
    y_min = 0.5*(sqrt(sqrt(2.*(a_t+b_t))+2.*sqrt(a_t))-sqrt(sqrt(2.*(a_t+b_t))-2.*sqrt(a_t)));
    y_max = pow(a_t,1./4.);
    closed_to_closed(yu[0],y_min,y_max,N[18]);
    //yb
    y_min = 0.5*(sqrt(2.*sqrt(b_t)+sqrt(2.*(a_t+b_t)))-sqrt(2.*sqrt(b_t)-sqrt(2.*(a_t+b_t))));
    y_max = pow(b_t,1./4.);
    closed_to_closed(yu[1],y_min,y_max,N[18]);
    
    for (i=0; i<N[18]; i++) {
        u_cmp[0][i] = gamma_III(yu[0][i],a_t,b_t,t_ps,t_th,'a','+','+');
        u_cmp[1][i] = gamma_III(yu[1][i],a_t,b_t,t_ps,t_th,'b','+','+');
    }
    u_cmp[0][N[18]-1].im = eps;
    u_cmp[1][N[18]-1].im = eps;
}

void build_grid_etap_eta2pi(double **s, int *N, double *s_below_cut, double *s_cv_plus, double *s_cv_minus, double *t_below_cut, double *t_cv_plus, double *t_cv_minus, bool pi0pi0) {
    int i,n,m;
    
    n = 0;
    m = 0;
    for (i=0; i<N[0]; i++) {
        s_below_cut[i+n] = s[0][i];
        t_below_cut[i+m] = s[0][i];
    }
    
    n = N[0];
    m = N[0];
    for (i=0; i<N[1]; i++) {
        s_below_cut[i+n] = s[1][i];
        t_below_cut[i+m] = s[1][i];
    }
    
    if (pi0pi0==true) {
        n = 0;
        m += N[1];
        for (i=0; i<N[2]; i++) {
            s_cv_plus[i+n] = s[2][i];
            t_below_cut[i+m] = s[2][i];
        }
        n += N[2];
    }
    else {
        n += N[1];
        m += N[1];
        for (i=0; i<N[2]; i++) {
            s_below_cut[i+n] = s[2][i];
            t_below_cut[i+m] = s[2][i];
        }
        n = 0;
    }
    
    m += N[2];
    for (i=0; i<N[3]; i++) {
        s_cv_plus[i+n] = s[3][i];
        t_below_cut[i+m] = s[3][i];
    }
    
    n += N[3];
    m += N[3];
    for (i=0; i<N[4]; i++) {
        s_cv_plus[i+n] = s[4][i];
        t_below_cut[i+m] = s[4][i];
    }
    
    n += N[4];
    m += N[4];
    for (i=0; i<N[5]; i++) {
        s_cv_plus[i+n] = s[5][i];
        t_below_cut[i+m] = s[5][i];
    }
    
    n += N[5];
    m = 0;
    for (i=0; i<N[6]; i++) {
        s_cv_plus[i+n] = s[6][i];
        t_cv_plus[i+m] = s[6][i];
    }
    
    n += N[6];
    m += N[6];
    for (i=0; i<N[7]; i++) {
        s_cv_plus[i+n] = s[7][i];
        t_cv_plus[i+m] = s[7][i];
    }
    
    n += N[7];
    m += N[7];
    for (i=0; i<N[8]; i++) {
        s_cv_plus[i+n] = s[8][i];
        t_cv_plus[i+m] = s[8][i];
    }
    
    n += N[8];
    m += N[8];
    for (i=0; i<N[9]; i++) {
        s_cv_plus[i+n] = s[9][i];
        t_cv_plus[i+m] = s[9][i];
    }
    
    n += N[9];
    m += N[9];
    for (i=0; i<N[10]; i++) {
        s_cv_plus[i+n] = s[10][i];
        t_cv_plus[i+m] = s[10][i];
    }
    
    n += N[10];
    m += N[10];
    for (i=0; i<N[11]; i++) {
        s_cv_plus[i+n] = s[11][i];
        t_cv_plus[i+m] = s[11][i];
    }
    
    n += N[11];
    m += N[11];
    for (i=0; i<N[12]; i++) {
        s_cv_plus[i+n] = s[12][i];
        t_cv_plus[i+m] = s[12][i];
    }
    
    n += N[12];
    m += N[12];
    for (i=0; i<N[13]; i++) {
        s_cv_plus[i+n] = s[13][i];
        t_cv_plus[i+m] = s[13][i];
    }
    
    for (i=0; i<N[14]; i++) {
        s_cv_minus[i] = s[14][i];
    }
    for (i=0; i<N[15]; i++) {
        t_cv_minus[i] = s[15][i];
    }
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

void build_s_complex(complex **s, double **y, double a, double b, double c, double d, char plusminus, char kappa_pm, int *N) {
    int i,n;
    double ya_min,ya_max,yb_min,yb_max,dya,dyb,eps;
    
    eps = 1.e-10; //for correct cv at the rim of the cut
    
    n = N[8];
    
    ya_min = 0.5*(sqrt(sqrt(2.*(a+b))+2.*sqrt(a))-sqrt(sqrt(2.*(a+b))-2.*sqrt(a)));
    ya_max = pow(a,1./4.);
    dya = (ya_max-ya_min)/(n-1);
    
    yb_min = 0.5*(sqrt(2.*sqrt(b)+sqrt(2.*(a+b)))-sqrt(2.*sqrt(b)-sqrt(2.*(a+b))));
    yb_max = pow(b,1./4.);
    dyb = (yb_max-yb_min)/(n-1);
    
    y[0][0] = ya_min;
    s[0][0] = gamma_III(y[0][0],a,b,c,d,'a',plusminus,kappa_pm);
    
    y[1][0] = yb_min;
    s[1][0] = gamma_III(y[1][0],a,b,c,d,'b',plusminus,kappa_pm);
    
    for (i=1; i<n-1; i++) {
        y[0][i] = ya_min+dya*i;
        s[0][i] = gamma_III(y[0][i],a,b,c,d,'a',plusminus,kappa_pm);
        
        y[1][i] = yb_min+dyb*i;
        s[1][i] = gamma_III(y[1][i],a,b,c,d,'b',plusminus,kappa_pm);
    }
    
    y[0][i] = ya_max;
    s[0][i] = gamma_III(y[0][i],a,b,c,d,'a',plusminus,kappa_pm);
    s[0][i].im = eps;
    
    y[1][i] = yb_max;
    s[1][i] = gamma_III(y[1][i],a,b,c,d,'b',plusminus,kappa_pm);
    s[1][i].im = eps;
}

void build_grid(double *s_cv_plus, double *s_cv_minus, double *s_below_cut, complex **s_cmp, double **y, double s_min, double s_max, double s_th, double eps_a, double eps_b, double a, double b, double c, double d, char plusminus, char kappa_pm, int *N, int Nsin) {
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
    
    build_s_complex(s_cmp,y,a,b,c,d,plusminus,kappa_pm,N);
    n = N[8];
    for (i=0; i<n; i++) {
        if (y[0][i]>=y[0][i+1] && i+1<n-1) {
            printf("FATAL ERROR: ya is not stictly monotonically increasing! step=%d, ya=%.10e\n\n",i,y[0][i]);
            exit(1);
        }
        if (y[1][i]>=y[1][i+1] && i+1<n-1) {
            printf("FATAL ERROR: yb is not stictly monotonically increasing! step=%d, yb=%.10e\n\n",i,y[1][i]);
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

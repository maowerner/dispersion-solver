#include "Basic.h"
#include "PhaseShifts.h"
#include "Grid.h"
#include "Omnes.h"
#include "AngularAverages.h"
#include "InhomogenitySqrt.h"
#include "InhomogenitySqrtCubed.h"
#include "Amplitude.h"
//#include "Iteration.h"
#include "InputOutput.h"

double test_func(double s) {
    return delta0(s)+delta_etapi(s);
}

void build_M_avg_u(complex_spline *M, complex *M_avg, double *s, int *N, double a, double b, double c, double d, double s0, char plusminus, char kappa_pm, int n) {
    M_avg_1(M[0],M_avg,s,a,b,c,d,plusminus,kappa_pm,N[0],n);
    M_avg_2(M[0],M[1],M_avg+N[0],s+N[0],a,b,c,d,s0,plusminus,kappa_pm,N[1],n);
    M_avg_III(M+7,M_avg+N[0]+N[1],s+N[0]+N[1],a,b,c,d,plusminus,kappa_pm,N[2],n);
    M_avg_4(M[2],M_avg+N[0]+N[1]+N[2],s+N[0]+N[1]+N[2],a,b,c,d,plusminus,kappa_pm,N[3],n);
}

double convergence_check(complex *Mi, complex *Mj, int N) {
    int i;
    double sum = 0.;
    
    for (i=0; i<N; i++) {
        sum += sqrt((pow(Mi[i].re-Mj[i].re,2.)+pow(Mi[i].im-Mj[i].im,2.)+1.e-10)/(Mi[i].re*Mi[i].re+Mi[i].im*Mi[i].im+1.e-10));
    }
    
    return sum/(double)N;
}

void copy_M(complex *Mi, complex *Mj, int N) {
    int i;
    
    for (i=0; i<N; i++) {
        Mj[i].re = Mi[i].re;
        Mj[i].im = Mi[i].im;
    }
}

int main (int argn, char **argc){
    int i,j,I,N[19],n,N_sin,n0,n1,n2,n0_sub,n1_sub,n2_sub;
    int N_amp_s[4],N_avg_s[4],N_amp_t[4],N_avg_t[4],N_amp_u,N_avg_u[4],N_hat_s[2],N_hat_t[2];
    double t,eps,eps_b,s_step,s_step_low_high,s,s0,L2,s_min,s_max,s_minus,t_minus,k,sum;
    double m1,m2,m3,m4; //particle masses
    double a,b,c,d;
    double a_s,b_s,s_ps,s_th,s_m,s_12,a_t,b_t,t_ps,t_th,t_m,t_12,a_u,b_u,u_ps,u_th,u_m,u_12;
    double *s_cv_plus,*s_cv_minus,*s_below_cut,*phi,*ya,*yb; //interpolation grid s-channel
    double *t_cv_plus,*t_cv_minus,*t_below_cut; //interpolation grid t-channel
    complex *s_complex,*sa_complex,*sb_complex; //interpolation grid
    complex F0[4],G0[4],F1[5],G1[5],F2[4],G2[4],z,Q12,Q32,R12,R32,temp0,temp1,temp2;
    char filename[1000];
    FILE *fin,*fout;
    
    gsl_set_error_handler_off();
    
    t = clock();
    
    //set particle masses
    
    //reading input file
    printf("Reading input file...\n");
    input(argc[1],&s_step,&L2,&n0_sub,&n1_sub,&n2_sub);
    printf("%.3f %d %d %d\n",s_step,n0_sub,n1_sub,n2_sub);
    printf("Done.\n\n");
    
    //pi0pi0 = false; //switch between pi0pi0 (true) and pi+pi- system (false)
    
    if (pi0pi0==true) {
        m1 = MPION0/MPION;
        m2 = MPION0/MPION;
    }
    else {
        m1 = 1.;
        m2 = 1.;
    }
    m3 = META/MPION;
    m4 = METAP/MPION;
    
    //s-channel variables
    s_th = pow(m1+m2,2.);
    s_ps = pow(m1-m2,2.);
    a_s = pow(m3-m4,2.);
    b_s = pow(m3+m4,2.);
    s_m = cmp_re(gamma_path(a_s,a_s,b_s,s_ps,s_th,'0','+'));
    s_12 = m1*(m4*m4-m3*m3)/(m3+m1);
    if (s_th>a_s) {
        printf("FATAL ERROR: particle masses have been set wrongly, m1<=m2<=m3<m4!\n");
        exit(1);
    }
    //printf("%.3f %.3f %.3f %.3f %.3f\n",s_th,a_s,b_s,s_m,s_12);
    
    //t-channel variables
    t_th = pow(m1+m3,2.);
    t_ps = pow(m1-m3,2.);
    a_t = pow(m2-m4,2.);
    b_t = pow(m2+m4,2.);
    t_m = cmp_re(gamma_path(a_t,a_t,b_t,t_ps,t_th,'-','-'));
    t_12 = 0.5*(m4*m4+m3*m3-2.*m1*m1);
    if (t_th>a_t) {
        printf("FATAL ERROR: particle masses have been set wrongly, m1<=m2<=m3<m4!\n");
        exit(1);
    }
    //printf("%.3f %.3f %.3f %.3f %.3f\n",t_th,a_t,b_t,t_m,t_12);
    
    //u-channel variables
    u_th = pow(m2+m3,2.);
    u_ps = pow(m2-m3,2.);
    a_u = pow(m1-m4,2.);
    b_u = pow(m1+m4,2.);
    u_m = cmp_re(gamma_path(a_u,a_u,b_u,u_ps,u_th,'+','+'));
    u_12 = (m3*(m4*m4-m1*m1)-m1*(m3*m3-m1*m1))/(m3+m1);
    if (u_th>a_u) {
        printf("FATAL ERROR: particle masses have been set wrongly, m1<=m2<=m3<m4!\n");
        exit(1);
    }
    //printf("%.3f %.3f %.3f %.3f %.3f\n\n",u_th,a_u,b_u,u_m,u_12);
    
    printf("Mass of decaying meson: %.3f MeV\n\n",METAP);
    
    //allocate gsl spline accelerators for the phase shifts
    delta0_params.acc = gsl_interp_accel_alloc();
    delta1_params.acc = gsl_interp_accel_alloc();
    delta2_params.acc = gsl_interp_accel_alloc();
    delta_etapi_params.acc = gsl_interp_accel_alloc();
    
    //import scattering phase shifts
    printf("Import scattering phase shifts...\n");
    sprintf(filename,"Input/Phases/%s",file_delta0);
    phase_input(filename,&delta0_params.spline);
    delta0_params.s_th = s_th;
    delta0_params = delta_high_energy_tale(delta0_params);
    sprintf(filename,"Input/Phases/%s",file_delta1);
    phase_input(filename,&delta1_params.spline);
    delta1_params.s_th = s_th;
    delta1_params = delta_high_energy_tale(delta1_params);
    sprintf(filename,"Input/Phases/%s",file_delta2);
    phase_input(filename,&delta2_params.spline);
    delta2_params.s_th = s_th;
    delta2_params = delta_high_energy_tale(delta2_params);
    sprintf(filename,"Input/Phases/%s",file_delta_etapi);
    phase_input(filename,&delta_etapi_params.spline);
    delta_etapi_params.s_th = t_th;
    delta_etapi_params = delta_high_energy_tale(delta_etapi_params);
    printf("Done.\n\n");
    
    //printf("a = %.3e\n",delta_etapi_params.a);
    //printf("b = %.3e\n",delta_etapi_params.b);
    //printf("L = %.3e\n",delta_etapi_params.L2);
    
    //s0 = s_pipi;
    //L2 = 1000.;
    
    printf("s_th = %.6f MPION**2\n",s_th);
    printf("t_th = %.6f MPION**2\n",t_th);
    printf("L2 = %.6f MPION**2\n\n",L2);
    
    s_max = L2;
    
    s_min = -s_max;
    
    eps = 1.0e-06;
    eps_b = 1.e-1; //treatment of dividing by 0 at b=METAP_P_MPION_SQUARED
    
    //for singularity in M_hat at a=(METAP-MPION)^2
    N_sin = (int)floor(log2(s_step/eps));
    //N_sin = 0;
    printf("N_singularity=%d %.2e\n\n",N_sin,s_step*pow(0.5,N_sin));
    
    s_step_low_high = 1.; //s_step for low and high energy ends
    //calculating number of grid points
    
    Ngrid_etap_eta2pi(N,N_sin,s_step,s_step_low_high,s_min,s_max,s_th,a_s,b_s,s_ps,s_m,s_12,t_th,a_t,b_t,t_ps,t_m,t_12,u_m,u_12,eps_b);
    
    //number of grid points for angular averages on the upper rim of the cut
    if (pi0pi0==true) {
        N_avg_s[0] = N[2]+N[3]; //[s_th,s_12]
        N_avg_s[1] = N[4]; //[s_12,a_s]
    }
    else {
        N_avg_s[0] = N[3]; //[s_th,s_12]
        N_avg_s[1] = N[4]; //[s_12,a_s]
    }
    
    N_avg_s[2] = N[5]+N[6]+N[7]+N[8]+N[9]+N[10]+N[11]+N[12]; //[a_s,b_s]
    N_avg_s[3] = N[13]; //[b_s,s_max]
    
    N_avg_t[0] = N[6]; //[t_th,t_12]
    N_avg_t[1] = N[7]+N[8]; //[t_12,a_t]
    N_avg_t[2] = N[9]+N[10]+N[11]; //[a_t,b_t]
    N_avg_t[3] = N[12]+N[13]; //[b_t,s_max]
    
    N_avg_u[0] = N[6]+N[7]; //[t_th,u_12]
    N_avg_u[1] = N[8]; //[t_12,a_t]
    N_avg_u[2] = N[9]+N[10]+N[11]; //[a_t,b_t]
    N_avg_u[3] = N[12]+N[13]; //[b_t,s_max]
    
    //number of grid points for Omnes-functions, amplitudes and inhomogeneities
    N_amp_s[0] = N_avg_s[0]+N_avg_s[1]+N_avg_s[2]+N_avg_s[3]; //upper rim of the cut [s_th,s_max]
    N_amp_s[1] = N[14]; //lower rim of the cut [s_th,t_m]
    if (pi0pi0==true) {
        N_amp_s[2] = N[0]+N[1]; //left from the cut [s_min,s_th]
    }
    else {
        N_amp_s[2] = N[0]+N[1]+N[2]; //left from the cut [s_min,s_th]
    }
    N_amp_s[3] = N[16]; //complex plane
    
    N_amp_t[0] = N_avg_t[0]+N_avg_t[1]+N_avg_t[2]+N_avg_t[3]; //upper rim of the cut [t_th,s_max]
    N_amp_t[1] = N[15]; //lower rim of the cut [t_th,s_m]
    N_amp_t[2] = N[0]+N[1]+N[2]+N[3]+N[4]+N[5]; //left from the cut [s_min,t_th]
    N_amp_t[3] = N[17]; //complex plane
    
    N_amp_u = N[18]; //complex plane
    
    //number of grid points for tilde and hat functions
    N_hat_s[0] = N_avg_s[0]+N_avg_s[1]; //[s_th,a_s]
    N_hat_s[1] = N_avg_s[2]+N_avg_s[3]; //[a_s,s_max]
    
    N_hat_t[0] = N_avg_t[0]+N_avg_t[1]; //[t_th,a_t]
    N_hat_t[1] = N_avg_t[2]+N_avg_t[3]; //[a_t,s_max]
    
    s_below_cut = (double*)malloc(N_amp_s[2]*sizeof(double));
    s_cv_plus = (double*)malloc(N_amp_s[0]*sizeof(double));
    s_cv_minus = (double*)malloc(N_amp_s[1]*sizeof(double));
    
    t_below_cut = (double*)malloc(N_amp_t[2]*sizeof(double));
    t_cv_plus = (double*)malloc(N_amp_t[0]*sizeof(double));
    t_cv_minus = (double*)malloc(N_amp_t[1]*sizeof(double));
    
    double **s_grid = (double**)malloc(16*sizeof(double*));
    for (i=0; i<16; i++) {
        s_grid[i] = (double*)malloc(N[i]*sizeof(double*));
    }
    
    double **ys = (double**)malloc(2*sizeof(double*));
    for (i=0; i<2; i++) {
        ys[i] = (double*)malloc(N_amp_s[3]*sizeof(double));
    }
    
    complex **s_cmp = (complex**)malloc(2*sizeof(complex*));
    for (i=0; i<2; i++) {
        s_cmp[i] = (complex*)malloc(N_amp_s[3]*sizeof(complex));
    }
    
    double **yt = (double**)malloc(2*sizeof(double*));
    for (i=0; i<2; i++) {
        yt[i] = (double*)malloc(N_amp_t[3]*sizeof(double));
    }
    
    complex **t_cmp = (complex**)malloc(2*sizeof(complex*));
    for (i=0; i<2; i++) {
        t_cmp[i] = (complex*)malloc(N_amp_t[3]*sizeof(complex));
    }
    
    double **yu = (double**)malloc(2*sizeof(double*));
    for (i=0; i<2; i++) {
        yu[i] = (double*)malloc(N_amp_u*sizeof(double));
    }
    
    complex **u_cmp = (complex**)malloc(2*sizeof(complex*));
    for (i=0; i<2; i++) {
        u_cmp[i] = (complex*)malloc(N_amp_u*sizeof(complex));
    }
    
    grid_etap_eta2pi(s_grid,ys,yt,yu,s_cmp,t_cmp,u_cmp,N,N_sin,s_step,s_step_low_high,s_min,s_max,s_th,a_s,b_s,s_ps,s_m,s_12,t_th,a_t,b_t,t_ps,t_m,t_12,u_m,u_12,eps_b,eps);
    
    build_grid_etap_eta2pi(s_grid,N,s_below_cut,s_cv_plus,s_cv_minus,t_below_cut,t_cv_plus,t_cv_minus,pi0pi0);
    
    //allocate memory for omnes_I functions
    //S-wave pi-pi
    complex **omnes0 = (complex **)malloc(5*sizeof(complex *));
    for (i=0; i<4; i++) {
        omnes0[i] = (complex *)malloc(N_amp_s[i]*sizeof(complex));
    }
    omnes0[4] = (complex *)malloc(N_amp_s[3]*sizeof(complex));
    //S-wave eta-pi
    complex **omnes1 = (complex **)malloc(7*sizeof(complex *));
    for (i=0; i<4; i++) {
        omnes1[i] = (complex *)malloc(N_amp_t[i]*sizeof(complex));
    }
    omnes1[4] = (complex *)malloc(N_amp_t[3]*sizeof(complex));
    omnes1[5] = (complex *)malloc(N_amp_u*sizeof(complex));
    omnes1[6] = (complex *)malloc(N_amp_u*sizeof(complex));
    
    
    //allocate memory for the inhomogeneity integrals
    complex **M0_inhom = (complex **)malloc(7*sizeof(complex *));
    for (i=0; i<4; i++) {
        M0_inhom[i] = (complex *)malloc(N_amp_s[i]*sizeof(complex));
    }
    for (i=4; i<7; i++) {
        M0_inhom[i] = (complex *)malloc(N_amp_s[3]*sizeof(complex));
    }
    
    complex **M1_inhom = (complex **)malloc(11*sizeof(complex *));
    for (i=0; i<4; i++) {
        M1_inhom[i] = (complex *)malloc(N_amp_t[i]*sizeof(complex));
    }
    for (i=4; i<7; i++) {
        M1_inhom[i] = (complex *)malloc(N_amp_t[3]*sizeof(complex));
    }
    for (i=7; i<11; i++) {
        M1_inhom[i] = (complex *)malloc(N_amp_u*sizeof(complex));
    }
    
    //allocate memory for the amplitudes
    complex_spline *M0 = (complex_spline *)malloc(7*sizeof(complex_spline));
    for (i=0; i<4; i++) {
        M0[i].re = gsl_spline_alloc(gsl_interp_cspline,N_amp_s[i]);
        M0[i].im = gsl_spline_alloc(gsl_interp_cspline,N_amp_s[i]);
    }
    for (i=4; i<7; i++) {
        M0[i].re = gsl_spline_alloc(gsl_interp_cspline,N_amp_s[3]);
        M0[i].im = gsl_spline_alloc(gsl_interp_cspline,N_amp_s[3]);
    }
    
    complex_spline *M1 = (complex_spline *)malloc(11*sizeof(complex_spline));
    for (i=0; i<4; i++) {
        M1[i].re = gsl_spline_alloc(gsl_interp_cspline,N_amp_t[i]);
        M1[i].im = gsl_spline_alloc(gsl_interp_cspline,N_amp_t[i]);
    }
    for (i=4; i<7; i++) {
        M1[i].re = gsl_spline_alloc(gsl_interp_cspline,N_amp_t[3]);
        M1[i].im = gsl_spline_alloc(gsl_interp_cspline,N_amp_t[3]);
    }
    for (i=7; i<11; i++) {
        M1[i].re = gsl_spline_alloc(gsl_interp_cspline,N_amp_u);
        M1[i].im = gsl_spline_alloc(gsl_interp_cspline,N_amp_u);
    }
    
    //allocate memory for the angular averages
    complex *M0_avg = (complex *)malloc(N_amp_s[0]*sizeof(complex));
    complex *M1_avg_s = (complex *)malloc(N_amp_t[0]*sizeof(complex));
    complex *M1_avg_u = (complex *)malloc(N_amp_t[0]*sizeof(complex));
    
    //allocate memory for the angular averages convergence check
    complex *M0_avg_convergence = (complex *)malloc(N_amp_s[0]*sizeof(complex));
    complex *M1_avg_s_convergence = (complex *)malloc(N_amp_t[0]*sizeof(complex));
    complex *M1_avg_u_convergence = (complex *)malloc(N_amp_t[0]*sizeof(complex));
    
    //allocate memory for the hat functions
    complex_spline *M0_hat = (complex_spline *)malloc(2*sizeof(complex_spline));
    complex_spline *M1_hat = (complex_spline *)malloc(2*sizeof(complex_spline));
    
    complex_spline_alloc(M0_hat,2,N_hat_s);
    complex_spline_alloc(M1_hat,2,N_hat_t);
    
    /*complex s_test,O_test;
    s_test.im = 1.e-10;
    for (i=0; i<8000; i++) {
        s_test.re = 0.1*i;
        O_test = Omnes_function(delta1,delta1_params,s_test);
        printf("%.6e %.6e %.6e %.6e %.6e\n",s_test.re*MPION*MPION/1.e+06,sqrt(O_test.re*O_test.re+O_test.im*O_test.im),delta1(s_test.re),O_test.re,O_test.im);
        //printf("%.6e %.6e\n",s_test.re,delta1(s_test.re));
    }
    exit(1);*/
    
    /*for (i=0; i<N_amp_s[2]; i++) {
        if (s_below_cut[i]>4.*pow(MPION0/MPION,2.)) {
            printf("%.10e %.10e\n",s_below_cut[i],0.,0.);
        }
    }
    for (i=0; i<N_amp_s[0]; i++) {
        printf("%.10e %.10e\n",s_cv_plus[i],delta0(s_cv_plus[i]));
    }
    //printf("%.6e\n",4.*pow(MPION0/MPION,2.));
    exit(1);*/
    
    printf("Computing Omnes functions...\n");
    printf("Isospin 0...\n");
    build_omnes(omnes0,delta0,delta0_params,s_cv_plus,s_cv_minus,s_below_cut,s_cmp,N_amp_s);
    printf("Done.\n");
    printf("Isospin 1...\n");
    build_omnes(omnes1,delta_etapi,delta_etapi_params,t_cv_plus,t_cv_minus,t_below_cut,t_cmp,N_amp_t);
    omnes_complex(omnes1[5],delta_etapi,delta_etapi_params,u_cmp[0],N_amp_u);
    omnes_complex(omnes1[6],delta_etapi,delta_etapi_params,u_cmp[1],N_amp_u);
    printf("Done.\n\n");
    
    double dO0,dO1,d2O0,d2O1;
    
    dO0 = dOmnes0(delta0,delta0_params);
    d2O0 = d2Omnes0(delta0,delta0_params,dO0);
    
    dO1 = dOmnes0(delta_etapi,delta_etapi_params);
    d2O1 = d2Omnes0(delta_etapi,delta_etapi_params,dO1);
    
    dO0 /= pow(MPION/1000,2.);
    d2O0 /= pow(MPION/1000,4.);
    
    dO1 /= pow(MPION/1000,2.);
    d2O1 /= pow(MPION/1000,4.);
    
    printf("Omnes0'(0) = %.10e, Omnes0''(0) = %.10e\n",dO0,d2O0);
    printf("Omnes1'(0) = %.10e, Omnes1''(0) = %.10e\n",dO1,d2O1);
    
    //printf("\ns_th = %.6e\n\n",delta0_params.s_th);
    
    //exit(1);
    
    /*for (i=0; i<N_amp_t[3]; i++) {
        printf("%.10e %.10e\n",t_cmp[0][i].re,t_cmp[0][i].im);
    }
    for (i=0; i<N_amp_t[3]; i++) {
        printf("%.10e %.10e\n",t_cmp[1][i].re,t_cmp[1][i].im);
    }
    exit(1);*/
    
    printf("Writing output for Omnes functions...\n");
    
    fout = fopen("OmnesFunctions","w");
    for (i=0; i<N_amp_s[2]; i++) {
        fprintf(fout,"%.10e %.10e %.10e %.10e %.10e\n",s_below_cut[i],omnes0[2][i].re,omnes0[2][i].im,omnes1[2][i].re,omnes1[2][i].im);
    }
    j = 0;
    for (i=0; i<N_amp_s[0]; i++) {
        if (s_cv_plus[i]<t_th) {
            fprintf(fout,"%.10e %.10e %.10e %.10e %.10e\n",s_cv_plus[i],omnes0[0][i].re,omnes0[0][i].im,omnes1[2][i+N_amp_s[2]].re,omnes1[2][i+N_amp_s[2]].im);
        }
        else {
            fprintf(fout,"%.10e %.10e %.10e %.10e %.10e\n",s_cv_plus[i],omnes0[0][i].re,omnes0[0][i].im,omnes1[0][j].re,omnes1[0][j].im);
            j++;
        }
    }
    /*for (i=0; i<N_amp_s[2]; i++) {
        fprintf(fout,"%.10e %.10e %.10e\n",s_below_cut[i],omnes0[2][i].re,omnes0[2][i].im);
    }
    j = 0;
    for (i=0; i<N_amp_s[0]; i++) {
        if (s_cv_plus[i]<t_th) {
            fprintf(fout,"%.10e %.10e %.10e\n",s_cv_plus[i],omnes0[0][i].re,omnes0[0][i].im);
        }
        else {
            fprintf(fout,"%.10e %.10e %.10e\n",s_cv_plus[i],omnes0[0][i].re,omnes0[0][i].im);
            j++;
        }
    }*/
    
    printf("Done.\n\n");
    exit(1);
    
    for (I=0; I<N_SUB_CONST; I++) {
        subtraction_constant(SUB_CONST[I],&n0,&n1,&n2);
        printf("Dertermine the system for subtraction constant '%s'\n\n",SUB_CONST[I]);
        
        printf("Setting initial inhomogeneities to zero...\n");
        for (i=0; i<3; i++) {
            for (j=0; j<N_amp_s[i]; j++) {
                M0_inhom[i][j].re = 0.;
                M0_inhom[i][j].im = 0.;
            }
        }
        for (i=3; i<7; i++) {
            for (j=0; j<N_amp_s[3]; j++) {
                M0_inhom[i][j].re = 0.;
                M0_inhom[i][j].im = 0.;
            }
        }
        
        for (i=0; i<3; i++) {
            for (j=0; j<N_amp_t[i]; j++) {
                M1_inhom[i][j].re = 0.;
                M1_inhom[i][j].im = 0.;
            }
        }
        for (i=3; i<7; i++) {
            for (j=0; j<N_amp_t[3]; j++) {
                M1_inhom[i][j].re = 0.;
                M1_inhom[i][j].im = 0.;
            }
        }
        for (i=7; i<11; i++) {
            for (j=0; j<N_amp_u; j++) {
                M1_inhom[i][j].re = 0.;
                M1_inhom[i][j].im = 0.;
            }
        }
        printf("Done.\n\n");
        
        printf("Computing Amplitudes...\n");
        build_amplitude(M0,omnes0,M0_inhom,s_cv_plus,s_cv_minus,s_below_cut,s_cmp,ys,s_th,L2,N_amp_s,n0);
        build_amplitude(M1,omnes1,M1_inhom,t_cv_plus,t_cv_minus,t_below_cut,t_cmp,yt,t_th,L2,N_amp_t,n1);
        amplitude_complex(M1[7],omnes1[5],M1_inhom[7],u_cmp[0],yu[0],'-',N_amp_u,n1);
        amplitude_complex(M1[8],omnes1[6],M1_inhom[8],u_cmp[1],yu[1],'-',N_amp_u,n1);
        amplitude_complex(M1[9],omnes1[5],M1_inhom[9],u_cmp[0],yu[0],'+',N_amp_u,n1);
        amplitude_complex(M1[10],omnes1[6],M1_inhom[10],u_cmp[1],yu[1],'+',N_amp_u,n1);
        printf("Done.\n\n");
        
        printf("Computing angular averages...\n");
        build_M_avg(M1,M0_avg,s_cv_plus,N_avg_s,a_s,b_s,s_ps,s_th,t_th,'0','+',0);
        build_M_avg(M0,M1_avg_s,t_cv_plus,N_avg_t,a_t,b_t,t_ps,t_th,s_th,'-','-',0);
        build_M_avg_u(M1,M1_avg_u,t_cv_plus,N_avg_u,a_t,b_t,t_ps,t_th,t_th,'+','+',0);
        printf("Done.\n\n");
        
        s = 50.0;
        sprintf(filename,"convergence_%s",SUB_CONST[I]);
        fout = fopen(filename,"w");
        gsl_interp_accel *acc0 = gsl_interp_accel_alloc();
        gsl_interp_accel *acc1 = gsl_interp_accel_alloc();
        
        sum = (convergence_check(M0_avg,M0_avg_convergence,N_amp_s[0])+convergence_check(M1_avg_s,M1_avg_s_convergence,N_amp_t[0])+convergence_check(M1_avg_u,M1_avg_u_convergence,N_amp_t[0]))/3.;
        
        copy_M(M0_avg,M0_avg_convergence,N_amp_s[0]);
        copy_M(M1_avg_s,M1_avg_s_convergence,N_amp_t[0]);
        copy_M(M1_avg_u,M1_avg_u_convergence,N_amp_t[0]);
        
        fprintf(fout,"%d %+.5e %+.5e %+.5e %+.5e %+.5e\n",0,gsl_spline_eval(M0[0].re,s,acc0),gsl_spline_eval(M0[0].im,s,acc0),gsl_spline_eval(M1[0].re,s,acc1),gsl_spline_eval(M1[0].im,s,acc1),sum);
        
        for (i=1; i<N_ITER; i++) {
            printf("Iteration %d:\n",i);
            
            printf("Building inhomogenity integrands...\n");
            build_inhomogenity_integrand_etap_eta2pi(M0_hat,M1_hat,F0,G0,F1,G1,M0_avg,M1_avg_s,M1_avg_u,omnes0[0],omnes1[0],s_cv_plus,t_cv_plus,s_th,s_ps,a_s,b_s,t_th,t_ps,a_t,b_t,L2,N_avg_s,N_avg_t,N_sin,n0_sub,n1_sub);
            printf("Done.\n\n");
            
            printf("Computing inhomogenities...\n");
            build_inhomogenity_sqrt(M0_inhom,M0_hat,F0,G0,s_th,a_s,L2,s_cv_plus,s_cv_minus,s_below_cut,s_cmp,N_amp_s,n0_sub);
            build_inhomogenity_sqrt(M1_inhom,M1_hat,F1,G1,t_th,a_t,L2,t_cv_plus,t_cv_minus,t_below_cut,t_cmp,N_amp_t,n1_sub);
            inhomogenity_sqrt_complex(M1_inhom[7],M1_hat,F1,G1,u_cmp[0],'-',t_th,a_t,L2,N_amp_u,n1_sub);
            inhomogenity_sqrt_complex(M1_inhom[8],M1_hat,F1,G1,u_cmp[1],'-',t_th,a_t,L2,N_amp_u,n1_sub);
            inhomogenity_sqrt_complex(M1_inhom[9],M1_hat,F1,G1,u_cmp[0],'+',t_th,a_t,L2,N_amp_u,n1_sub);
            inhomogenity_sqrt_complex(M1_inhom[10],M1_hat,F1,G1,u_cmp[1],'+',t_th,a_t,L2,N_amp_u,n1_sub);
            printf("Done.\n\n");
            
            printf("Computing Amplitudes...\n");
            build_amplitude(M0,omnes0,M0_inhom,s_cv_plus,s_cv_minus,s_below_cut,s_cmp,ys,s_th,L2,N_amp_s,n0);
            build_amplitude(M1,omnes1,M1_inhom,t_cv_plus,t_cv_minus,t_below_cut,t_cmp,yt,t_th,L2,N_amp_t,n1);
            amplitude_complex(M1[7],omnes1[5],M1_inhom[7],u_cmp[0],yu[0],'-',N_amp_u,n1);
            amplitude_complex(M1[8],omnes1[6],M1_inhom[8],u_cmp[1],yu[1],'-',N_amp_u,n1);
            amplitude_complex(M1[9],omnes1[5],M1_inhom[9],u_cmp[0],yu[0],'+',N_amp_u,n1);
            amplitude_complex(M1[10],omnes1[6],M1_inhom[10],u_cmp[1],yu[1],'+',N_amp_u,n1);
            printf("Done.\n\n");
            
            printf("Computing angular averages...\n");
            build_M_avg(M1,M0_avg,s_cv_plus,N_avg_s,a_s,b_s,s_ps,s_th,t_th,'0','+',0);
            build_M_avg(M0,M1_avg_s,t_cv_plus,N_avg_t,a_t,b_t,t_ps,t_th,s_th,'-','-',0);
            build_M_avg_u(M1,M1_avg_u,t_cv_plus,N_avg_u,a_t,b_t,t_ps,t_th,t_th,'+','+',0);
            printf("Done.\n\n");
            
            sum = (convergence_check(M0_avg,M0_avg_convergence,N_amp_s[0])+convergence_check(M1_avg_s,M1_avg_s_convergence,N_amp_t[0])+convergence_check(M1_avg_u,M1_avg_u_convergence,N_amp_t[0]))/3.;
            
            copy_M(M0_avg,M0_avg_convergence,N_amp_s[0]);
            copy_M(M1_avg_s,M1_avg_s_convergence,N_amp_t[0]);
            copy_M(M1_avg_u,M1_avg_u_convergence,N_amp_t[0]);
            
            s = 50.0;
            fprintf(fout,"%d %+.5e %+.5e %+.5e %+.5e %+.5e\n",i,gsl_spline_eval(M0[0].re,s,acc0),gsl_spline_eval(M0[0].im,s,acc0),gsl_spline_eval(M1[0].re,s,acc1),gsl_spline_eval(M1[0].im,s,acc1),sum);
        }
        fclose(fout);
        gsl_interp_accel_free(acc0);
        gsl_interp_accel_free(acc1);
        
        printf("Writing output for subraction constant '%s'...\n",SUB_CONST[I]);
        
        sprintf(filename,"M_tilde_%s",SUB_CONST[I]);
        fout = fopen(filename,"w");
        j = 0;
        for (i=0; i<N_amp_s[0]; i++) {
            if (s_cv_plus[i]<t_th) {
                fprintf(fout,"%.10e %.10e %.10e\n",s_cv_plus[i],M0_avg[i].re,M0_avg[i].im);
            }
            else {
                fprintf(fout,"%.10e %.10e %.10e %.10e %.10e %.10e %.10e\n",s_cv_plus[i],M0_avg[i].re,M0_avg[i].im,M1_avg_s[j].re,M1_avg_s[j].im,M1_avg_u[j].re,M1_avg_u[j].im);
                j++;
            }
        }
        fclose(fout);
        
        sprintf(filename,"M_%s",SUB_CONST[I]);
        fout = fopen(filename,"w");
        gsl_interp_accel *acc_s0 = gsl_interp_accel_alloc();
        gsl_interp_accel *acc_t0 = gsl_interp_accel_alloc();
        gsl_interp_accel *acc_s1 = gsl_interp_accel_alloc();
        gsl_interp_accel *acc_t1 = gsl_interp_accel_alloc();
        for (i=0; i<N_amp_s[2]; i++) {
            fprintf(fout,"%.10e %.10e %.10e %.10e %.10e\n",s_below_cut[i],gsl_spline_eval(M0[2].re,s_below_cut[i],acc_s0),gsl_spline_eval(M0[2].im,s_below_cut[i],acc_s0),gsl_spline_eval(M1[2].re,s_below_cut[i],acc_t0),gsl_spline_eval(M1[2].im,s_below_cut[i],acc_t0));
        }
        for (i=0; i<N_amp_s[0]; i++) {
            if (s_cv_plus[i]<t_th) {
                fprintf(fout,"%.10e %.10e %.10e %.10e %.10e\n",s_cv_plus[i],gsl_spline_eval(M0[0].re,s_cv_plus[i],acc_s1),gsl_spline_eval(M0[0].im,s_cv_plus[i],acc_s1),gsl_spline_eval(M1[2].re,s_cv_plus[i],acc_t0),gsl_spline_eval(M1[2].im,s_cv_plus[i],acc_t0));
            }
            else {
                fprintf(fout,"%.10e %.10e %.10e %.10e %.10e\n",s_cv_plus[i],gsl_spline_eval(M0[0].re,s_cv_plus[i],acc_s1),gsl_spline_eval(M0[0].im,s_cv_plus[i],acc_s1),gsl_spline_eval(M1[0].re,s_cv_plus[i],acc_t1),gsl_spline_eval(M1[0].im,s_cv_plus[i],acc_t1));
            }
        }
        fclose(fout);
        
        sprintf(filename,"integrand_%s",SUB_CONST[I]);
        gsl_interp_accel *acc_s2 = gsl_interp_accel_alloc();
        gsl_interp_accel *acc_t2 = gsl_interp_accel_alloc();
        gsl_interp_accel *acc_s3 = gsl_interp_accel_alloc();
        gsl_interp_accel *acc_t3 = gsl_interp_accel_alloc();
        fout = fopen(filename,"w");
        for (i=0; i<N_amp_s[0]; i++) {
            if (s_cv_plus[i]<t_th && s_cv_plus[i]<a_s) {
                fprintf(fout,"%.10e %.10e %.10e %.10e %.10e\n",s_cv_plus[i],gsl_spline_eval(M0_hat[0].re,s_cv_plus[i],acc_s2),gsl_spline_eval(M0_hat[0].im,s_cv_plus[i],acc_s2),0.,0.);
            }
            else if (s_cv_plus[i]<t_th && s_cv_plus[i]>a_s) {
                fprintf(fout,"%.10e %.10e %.10e %.10e %.10e\n",s_cv_plus[i],gsl_spline_eval(M0_hat[1].re,s_cv_plus[i],acc_s3),gsl_spline_eval(M0_hat[1].im,s_cv_plus[i],acc_s3),0.,0.);
            }
            else if (s_cv_plus[i]>t_th && s_cv_plus[i]<a_t) {
                fprintf(fout,"%.10e %.10e %.10e %.10e %.10e\n",s_cv_plus[i],gsl_spline_eval(M0_hat[1].re,s_cv_plus[i],acc_s3),gsl_spline_eval(M0_hat[1].im,s_cv_plus[i],acc_s3),gsl_spline_eval(M1_hat[0].re,s_cv_plus[i],acc_t2),gsl_spline_eval(M1_hat[0].im,s_cv_plus[i],acc_t2));
            }
            else if (s_cv_plus[i]>t_th && s_cv_plus[i]>a_t) {
                fprintf(fout,"%.10e %.10e %.10e %.10e %.10e\n",s_cv_plus[i],gsl_spline_eval(M0_hat[1].re,s_cv_plus[i],acc_s3),gsl_spline_eval(M0_hat[1].im,s_cv_plus[i],acc_s3),gsl_spline_eval(M1_hat[1].re,s_cv_plus[i],acc_t3),gsl_spline_eval(M1_hat[1].im,s_cv_plus[i],acc_t3));
            }
        }
        fclose(fout);
        
        sprintf(filename,"M_int_%s",SUB_CONST[I]);
        fout = fopen(filename,"w");
        for (i=0; i<N_amp_s[2]; i++) {
            fprintf(fout,"%.10e %.10e %.10e %.10e %.10e\n",s_below_cut[i],M0_inhom[2][i].re,M0_inhom[2][i].im,M1_inhom[2][i].re,M1_inhom[2][i].im);
        }
        j = 0;
        for (i=0; i<N_amp_s[0]; i++) {
            if (s_cv_plus[i]<t_th) {
                fprintf(fout,"%.10e %.10e %.10e %.10e %.10e\n",s_cv_plus[i],M0_inhom[0][i].re,M0_inhom[0][i].im,M1_inhom[2][i+N_amp_s[2]].re,M1_inhom[2][i+N_amp_s[2]].im);
            }
            else {
                fprintf(fout,"%.10e %.10e %.10e %.10e %.10e\n",s_cv_plus[i],M0_inhom[0][i].re,M0_inhom[0][i].im,M1_inhom[0][j].re,M1_inhom[0][j].im);
                j++;
            }
        }
        
        complex test_tmp;
        
        test_tmp = inhomogenity_sqrt_below_cut_single_value(M0_hat,F0,G0,0.,s_th,a_s,b_s,n0_sub);
        printf("I0_%s = %+.6e %+.6e i\n",SUB_CONST[I],test_tmp.re,test_tmp.im);
        
        test_tmp = inhomogenity_sqrt_below_cut_single_value(M1_hat,F1,G1,0.,t_th,a_t,b_t,n0_sub);
        printf("I1_%s = %+.6e %+.6e i\n",SUB_CONST[I],test_tmp.re,test_tmp.im);
        
        printf("Done.\n\n");
    }
    
    
//
//    free(s_cv_plus);
//    free(s_cv_minus);
//    free(s_complex);
//    
//    for (i=0; i<4; i++) {
//        free(omnes0[i]);
//    }
//    free(omnes0);
//    
//    for (i=0; i<4; i++) {
//        free(omnes1[i]);
//    }
//    free(omnes1);
//    
//    for (i=0; i<4; i++) {
//        free(omnes2[i]);
//    }
//    free(omnes2);
//    
//    for (i=0; i<4; i++) {
//        free(M0_inhom[i]);
//    }
//    free(M0_inhom);
//    
//    for (i=0; i<4; i++) {
//        free(M1_inhom[i]);
//    }
//    free(M1_inhom);
//    
//    for (i=0; i<4; i++) {
//        free(M2_inhom[i]);
//    }
//    free(M2_inhom);
//    
//    n=2;
//    for (i=0; i<n; i++) {
//        free(M0_avg[i]);
//    }
//    free(M0_avg);
//    
//    n=3;
//    for (i=0; i<n; i++) {
//        free(M1_avg[i]);
//    }
//    free(M1_avg);
//    
//    n=2;
//    for (i=0; i<n; i++) {
//        free(M2_avg[i]);
//    }
//    free(M2_avg);
//    
//    complex_spline_free(M0_hat,2);
//    complex_spline_free(M1_hat,2);
//    complex_spline_free(M2_hat,2);
//    free(M0_hat);
//    free(M1_hat);
//    free(M2_hat);
//    
//    complex_spline_free(M0,4);
//    complex_spline_free(M1,4);
//    complex_spline_free(M2,4);
//    free(M0);
//    free(M1);
//    free(M2);
    
    gsl_spline_free(delta0_params.spline);
    gsl_spline_free(delta1_params.spline);
    gsl_spline_free(delta2_params.spline);
    gsl_spline_free(delta_etapi_params.spline);
    
    gsl_interp_accel_free(delta0_params.acc);
    gsl_interp_accel_free(delta1_params.acc);
    gsl_interp_accel_free(delta2_params.acc);
    gsl_interp_accel_free(delta_etapi_params.acc);
    
    t = (clock()-t)/CLOCKS_PER_SEC;
    
    printf("\ncalculation time = %.3fs\n\n", t);
    
    return 0;
}

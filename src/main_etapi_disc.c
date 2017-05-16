#include "Basic.h"
#include "Grid.h"
#include "Omnes.h"
#include "AngularAverages.h"
#include "Etapi_disc.h"
#include "InhomogenitySqrt.h"
#include "InhomogenitySqrtCubed.h"
#include "Amplitude.h"
#include "InputOutput.h"

int main (int argn, char **argc){
    int i,j,N[10],n,m,N_sin,n0,n1,n2,n0_sub,n1_sub,n2_sub;
    int N_M_int[4],N_M_avg[4],N_M_hat[2],N_temp[2];
    double t,si,sf,ds,eps,eps_b,s_step,s,s0,L2,s_min,s_max,k,s_step_low_high;
    double *s_cv_plus,*s_cv_minus,*s_below_cut,*phi,*s_etapi_disc; //interpolation grid
    complex *s_complex; //interpolation grid
    complex F0[4],G0[4],F1[5],G1[5],F2[4],G2[4],z,Q12,Q32,R12,R32,temp0,temp1,temp2;
    char sub_const[2],phase_type[100],filename[1000];
    FILE *fin,*fout;
    
    gsl_set_error_handler_off();
    
    t = clock();
    
    printf("Mass of decaying meson: %.3f MeV\n\n",METAP);
    
    //reading input file
    printf("Reading input file...\n");
    Input(argc[1],&s_step,sub_const,&n0_sub,&n1_sub,&n2_sub,phase_type);
    //printf("%.3f %s %d %d %d\n",s_step,sub_const,n0_sub,n1_sub,n2_sub);
    //printf("%s\n",filename);
    printf("Done.\n\n");
    
    printf("Integration-grid-spacing: %.3fmpi**2\n\n",s_step);
    
    subtraction_constant(sub_const,&n0,&n1,&n2);
    printf("Subtractions: n0=%d, n1=%d, n2=%d\n\n",n0_sub,n1_sub,n2_sub);
    printf("Dertermine the system for subtraction constant '%s'\n\n",sub_const);
    
    gsl_spline *delta0;
    gsl_spline *delta1;
    gsl_spline *delta2;
    
    //import phases
    printf("Import scattering phase shifts...\n");
    sprintf(filename,"Input/Phases/%s_phases",phase_type);
    Phases(filename,&delta0,&delta1,&delta2,&s0,&L2);
    printf("Done.\n\n");
    
    //s0 = s_etapi;
    printf("s_th = %.6f MPION**2\n",s0);
    printf("L2 = %.6f MPION**2\n\n",L2);
    
    //import partial waves
    printf("Import partial waves...\n");
    
    complex_spline f0etapi,f2etapi,Metapi,*Metapi_tilde;
    
    n = 19513;//filesize f_eta->3pi
    m = 47625;//filesize f_etap->eta2pi
    i = 535;//point of s_sing in f_etap->eta2pi
    
    f0etapi.re = gsl_spline_alloc(gsl_interp_cspline,n);
    f0etapi.im = gsl_spline_alloc(gsl_interp_cspline,n);
    f2etapi.re = gsl_spline_alloc(gsl_interp_cspline,n);
    f2etapi.im = gsl_spline_alloc(gsl_interp_cspline,n);
    Metapi.re = gsl_spline_alloc(gsl_interp_cspline,m);
    Metapi.im = gsl_spline_alloc(gsl_interp_cspline,m);
    
    Metapi_tilde = (complex_spline*)malloc(2*sizeof(complex_spline));
    Metapi_tilde[0].re = gsl_spline_alloc(gsl_interp_cspline,i);
    Metapi_tilde[0].im = gsl_spline_alloc(gsl_interp_cspline,i);
    Metapi_tilde[1].re = gsl_spline_alloc(gsl_interp_cspline,m-i);
    Metapi_tilde[1].im = gsl_spline_alloc(gsl_interp_cspline,m-i);
    
    sprintf(filename,"../etapi_disc/f0_etapi_pipi.input");
    Partial_wave_eta_3pi_input(filename,f0etapi,L2,n);
    sprintf(filename,"../etapi_disc/f2_etapi_pipi.input");
    Partial_wave_eta_3pi_input(filename,f2etapi,L2,n);
    sprintf(filename,"../etapi_disc/f0_etappi_etapi.input");
    Partial_wave_etap_eta2pi_input(filename,Metapi,Metapi_tilde,L2,m,i);
    sprintf(filename,"test_partialwave");
    Output_test(filename,f0etapi,f2etapi,Metapi,Metapi_tilde,0.05,L2);
    printf("Done.\n\n");
    
    printf("L2 = %.6f MPION**2\n\n",L2);
    
//    s_max = L2;
//    s_min = integration_s_plus(s_max);
//    
//    eps = 1.0e-06;//treatment of dividing by 0 at a=METAP_M_MPION_SQUARED
//    eps_b = 1.0e-06; //treatment of dividing by 0 at b=METAP_P_MPION_SQUARED
//    
//    //for singularity in M_hat at a=(METAP-MPION)^2
//    N_sin = (int)floor(log2(s_step/eps));
//    //N_sin = 0;
//    printf("N_singularity=%d %.2e\n\n",N_sin,s_step*pow(0.5,N_sin));
//    
//    s_step_low_high = 1.; //s_step for low and high energy ends
//    //calculating number of grid points
//    grid_point_size(N,N_sin,s_step,s_step_low_high,s_min,s_max,s0,eps_b);
//    
//    //number of grid points for angular averages
//    N_M_avg[0] = N[2];
//    N_M_avg[1] = N[3]+N[9];
//    N_M_avg[2] = N[4];
//    N_M_avg[3] = N[5]+N[6];
//    
//    //number of grid points for Omnes-functions, amplitudes and inhomogeneities
//    N_M_int[0] = N_M_avg[0]+N_M_avg[1]+N_M_avg[2]+N_M_avg[3];
//    N_M_int[1] = N[7];
//    N_M_int[2] = N[8];
//    N_M_int[3] = N[0]+N[1];
//    
//    N_M_hat[0] = (int)ceil((METAP_M_MPION_SQUARED-s_etapi)/s_step)+N_sin;
//    N_temp[0] = (int)ceil((100.-METAP_M_MPION_SQUARED)/s_step)+N_sin;
//    N_temp[1] = (int)ceil((s_max-100.)/s_step_low_high);
//    N_M_hat[1] = N_temp[0]+N_temp[1];
//    
//    printf("Number of grid points for the Omnes functions and amplitudes:\n");
//    for (i=0; i<4; i++) {
//        printf("N_M_int(%d)=%d\n",i+1,N_M_int[i]);
//    }
//    printf("\n");
//    
//    printf("Number of grid points for angular averages:\n");
//    for (i=0; i<4; i++) {
//        printf("N_M_avg(%d)=%d\n",i+1,N_M_avg[i]);
//    }
//    printf("\n");
//    
//    printf("Number of grid points for inhomogenities:\n");
//    for (i=0; i<2; i++) {
//        printf("N_M_hat(%d)=%d\n",i+1,N_M_hat[i]);
//    }
//    printf("\n");
//    
//    s_below_cut = (double*)malloc(N_M_int[3]*sizeof(double));
//    s_cv_plus = (double*)malloc(N_M_int[0]*sizeof(double));
//    s_cv_minus = (double*)malloc(N_M_int[1]*sizeof(double));
//    phi = (double*)malloc(N_M_int[2]*sizeof(double));
//    s_complex = (complex*)malloc(N_M_int[2]*sizeof(complex));
//    s_etapi_disc = (double *)malloc((N_M_hat[0]+N_M_hat[1])*sizeof(double));
//    
//    build_grid(s_cv_plus,s_cv_minus,s_below_cut,s_complex,phi,s_min,s_max,s0,eps,eps_b,N,N_sin);
//    build_s_etapi_disc(s_etapi_disc,s_etapi,L2,eps,N_M_hat,N_temp,N_sin);
//    
////    for (i=0; i<(N_M_hat[0]+N_M_hat[1]); i++) {
////        printf("%d %.6e\n",i,s_etapi_disc[i]);
////        if (s_etapi_disc[i]>=s_etapi_disc[i+1] && i+1<(N_M_hat[0]+N_M_hat[1])) {
////            printf("nan: %d %.6e\n",i,s_etapi_disc[i]);
////        }
////    }
////    exit(1);
//    
//    complex *R12_cv_plus = (complex *)malloc(N_M_int[0]*sizeof(complex));
//    complex *Q12_cv_plus = (complex *)malloc(N_M_int[0]*sizeof(complex));
//    
//    for (i=0; i<N_M_int[0]; i++) {
//        R12_cv_plus[i] = R_sqrt(s_cv_plus[i],s_etapi,L2);
//        
//        Q12_cv_plus[i] = Q_sqrt(s_cv_plus[i],s_etapi,L2);
//    }
//    
//    complex *R12_cv_minus = (complex *)malloc(N_M_int[1]*sizeof(complex));
//    complex *Q12_cv_minus = (complex *)malloc(N_M_int[1]*sizeof(complex));
//    
//    for (i=0; i<N_M_int[1]; i++) {
//        R12_cv_minus[i] = R_sqrt(s_cv_minus[i],s_etapi,L2);
//        
//        Q12_cv_minus[i] = Q_sqrt(s_cv_minus[i],s_etapi,L2);
//    }
//    
//    complex *Q12_complex_f = (complex *)malloc(N_M_int[2]*sizeof(complex));
//    complex *Q12_complex_g = (complex *)malloc(N_M_int[2]*sizeof(complex));
//    
//    for (i=0; i<N_M_int[2]; i++) {
//        Q12_complex_f[i] = Q_sqrt_complex_f(s_complex[i],s_etapi);
//        Q12_complex_g[i] = Q_sqrt_complex_g(s_complex[i],L2);
//    }
//    
//    complex *Q12_below_cut = (complex *)malloc(N_M_int[3]*sizeof(complex));
//    
//    for (i=0; i<N_M_int[3]; i++) {
//        Q12_below_cut[i] = Q_sqrt(s_below_cut[i],s_etapi,L2);
//    }
//    
//    double *abs_omnes0 = (double *)malloc((N_M_hat[0]+N_M_hat[1])*sizeof(double));
//    double *abs_omnes2 = (double *)malloc((N_M_hat[0]+N_M_hat[1])*sizeof(double));
//    
//    double *abs_omnes_dump = (double *)malloc(N_M_int[0]*sizeof(double));
//    
//    complex **omnes0 = (complex **)malloc(4*sizeof(complex *));
//    for (i=0; i<4; i++) {
//        omnes0[i] = (complex *)malloc(N_M_int[i]*sizeof(complex));
//    }
//    
//    complex **omnes1 = (complex **)malloc(4*sizeof(complex *));
//    for (i=0; i<4; i++) {
//        omnes1[i] = (complex *)malloc(N_M_int[i]*sizeof(complex));
//    }
//    
//    complex **omnes2 = (complex **)malloc(4*sizeof(complex *));
//    for (i=0; i<4; i++) {
//        omnes2[i] = (complex *)malloc(N_M_int[i]*sizeof(complex));
//    }
//    
//    n=2;
//    complex **M0_avg = (complex **)malloc(n*sizeof(complex *));
//    for (i=0; i<n; i++) {
//        M0_avg[i] = (complex *)malloc(N_M_int[0]*sizeof(complex));
//    }
//    
//    n=3;
//    complex **M1_avg = (complex **)malloc(n*sizeof(complex *));
//    for (i=0; i<n; i++) {
//        M1_avg[i] = (complex *)malloc(N_M_int[0]*sizeof(complex));
//    }
//    
//    n=2;
//    complex **M2_avg = (complex **)malloc(n*sizeof(complex *));
//    for (i=0; i<n; i++) {
//        M2_avg[i] = (complex *)malloc(N_M_int[0]*sizeof(complex));
//    }
//    
//    complex_spline M0_finite;
//    complex_spline M2_finite;
//    complex_spline *M0_sing = (complex_spline *)malloc(2*sizeof(complex_spline));
//    complex_spline *M2_sing = (complex_spline *)malloc(2*sizeof(complex_spline));
//    
//    M0_finite.re = gsl_spline_alloc(gsl_interp_cspline,(N_M_hat[0]+N_M_hat[1]));
//    M0_finite.im = gsl_spline_alloc(gsl_interp_cspline,(N_M_hat[0]+N_M_hat[1]));
//    M2_finite.re = gsl_spline_alloc(gsl_interp_cspline,(N_M_hat[0]+N_M_hat[1]));
//    M2_finite.im = gsl_spline_alloc(gsl_interp_cspline,(N_M_hat[0]+N_M_hat[1]));
//    complex_spline_alloc(M0_sing,2,N_M_hat);
//    complex_spline_alloc(M2_sing,2,N_M_hat);
//    
//    complex_spline *M0 = (complex_spline *)malloc(4*sizeof(complex_spline));
//    complex_spline *M1 = (complex_spline *)malloc(4*sizeof(complex_spline));
//    complex_spline *M2 = (complex_spline *)malloc(4*sizeof(complex_spline));
//    
//    complex_spline_alloc(M0,4,N_M_int);
//    complex_spline_alloc(M1,4,N_M_int);
//    complex_spline_alloc(M2,4,N_M_int);
//    
//    printf("Computing Omnes functions...\n");
//    abs_omnes_cv_plus(abs_omnes0,delta0,s0,L2,2.*M_PI,s_etapi_disc,(N_M_hat[0]+N_M_hat[1]));
//    abs_omnes_cv_plus(abs_omnes2,delta2,s0,L2,0.,s_etapi_disc,(N_M_hat[0]+N_M_hat[1]));
//    build_omnes(omnes0,abs_omnes_dump,delta0,s0,L2,2.*M_PI,s_cv_plus,s_cv_minus,s_below_cut,s_complex,N_M_int);
//    build_omnes(omnes1,abs_omnes_dump,delta1,s0,L2,M_PI,s_cv_plus,s_cv_minus,s_below_cut,s_complex,N_M_int);
//    build_omnes(omnes2,abs_omnes_dump,delta2,s0,L2,0.,s_cv_plus,s_cv_minus,s_below_cut,s_complex,N_M_int);
//    printf("Done.\n\n");
//    
//    printf("Building inhomogenity integrands...\n");
//    build_inhomogenity_integrand_etapi_disc(M0_finite,M2_finite,M0_sing,M2_sing,F0,G0,F2,G2,f0etapi,f2etapi,Metapi,Metapi_tilde,delta0,delta2,abs_omnes0,abs_omnes2,s_etapi_disc,s_etapi,L2,N_M_hat,N_sin,n0_sub,n2_sub);
//    printf("Done.\n\n");
//    
//    complex **M0_inhom = (complex **)malloc(4*sizeof(complex *));
//    for (i=0; i<4; i++) {
//        M0_inhom[i] = (complex *)malloc(N_M_int[i]*sizeof(complex));
//    }
//    
//    complex **M1_inhom = (complex **)malloc(4*sizeof(complex *));
//    for (i=0; i<4; i++) {
//        M1_inhom[i] = (complex *)malloc(N_M_int[i]*sizeof(complex));
//    }
//    
//    complex **M2_inhom = (complex **)malloc(4*sizeof(complex *));
//    for (i=0; i<4; i++) {
//        M2_inhom[i] = (complex *)malloc(N_M_int[i]*sizeof(complex));
//    }
//    
//    printf("Calculate etapi-inhomogeneities...\n");
//    build_etapi_inhomogeneity(M0_inhom,M1_inhom,M2_inhom,M0_finite,M2_finite,M0_sing,M2_sing,F0,F2,G0,G2,s_cv_plus,s_cv_minus,s_below_cut,s_complex,s_min,L2,N_M_int,n0_sub,n1_sub,n2_sub);
//    printf("Done.\n\n");
//    
//    printf("Computing Amplitudes...\n");
//    build_amplitude(M0,omnes0,M0_inhom,s_cv_plus,s_cv_minus,s_below_cut,s_complex,phi,s0,L2,s_min,N_M_int,-1);
//    build_amplitude(M1,omnes1,M1_inhom,s_cv_plus,s_cv_minus,s_below_cut,s_complex,phi,s0,L2,s_min,N_M_int,-1);
//    build_amplitude(M2,omnes2,M2_inhom,s_cv_plus,s_cv_minus,s_below_cut,s_complex,phi,s0,L2,s_min,N_M_int,-1);
//    printf("Done.\n\n");
//    
//    printf("Computing angular averages...\n");
//    for (n=0; n<2; n++) {
//        build_M_avg(M0,M0_avg[n],s_cv_plus,N_M_avg,n);
//        build_M_avg(M1,M1_avg[n],s_cv_plus,N_M_avg,n);
//        build_M_avg(M2,M2_avg[n],s_cv_plus,N_M_avg,n);
//    }
//    build_M_avg(M1,M1_avg[n],s_cv_plus,N_M_avg,n);
//    printf("Done.\n\n");
//    
//    printf("Writing output...\n");
//    
////    sprintf(filename,"F_%s","etapi");
////    Amplitudes_output(filename,M0,M1,M2,s_below_cut,s_cv_plus,N_M_int[3],N_M_int[0]);
//    
//    sprintf(filename,"F_%s",sub_const);
//    Dispersion_integral_output(filename,M0_inhom,M1_inhom,M2_inhom,s_below_cut,s_cv_plus,N_M_int[3],N_M_int[0]);
//    
//    sprintf(filename,"A_%s","etapi");
//    Inhomogenities_output(filename,M0_avg,M1_avg,M2_avg,s_cv_plus,N_M_int[0]);
    
    printf("Done.\n\n");
    
    t = (clock()-t)/CLOCKS_PER_SEC;
    
    printf("\ncalculation time = %.3fs\n\n", t);
    
    return 0;
}

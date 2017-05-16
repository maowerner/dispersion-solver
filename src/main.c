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

int main (int argn, char **argc){
    int i,j,N[10],n,N_sin,n0,n1,n2,n0_sub,n1_sub,n2_sub;
    int N_M_int[7],N_M_avg[4],N_M_hat[2];
    double t,si,sf,ds,eps,eps_b,s_step,s_step_low_high,s,s0,L2,s_min,s_max,k;
    double a,b,c,d; //particle masses: a=(m1-m2)**2, b=(m1+m2)**2, c=(m3-m4)**2, d=(m3+m4)**2
    double *s_cv_plus,*s_cv_minus,*s_below_cut,**ys,*ya,*yb; //interpolation grid
    complex *s_complex,*sa_complex,*sb_complex,**s_cmp; //interpolation grid
    complex F0[4],G0[4],F1[5],G1[5],F2[4],G2[4],z,Q12,Q32,R12,R32,temp0,temp1,temp2;
    char sub_const[2],**file_phases,filename[1000];
    FILE *fin,*fout;
    
    gsl_set_error_handler_off();
    
    t = clock();
    
    a = METAP_M_MPION_SQUARED;
    b = METAP_P_MPION_SQUARED;
    c = 0.;
    d = 4.;
    
    printf("Mass of decaying meson: %.3f MeV\n\n",METAP);
    
    file_phases = (char**)malloc(4*sizeof(char*));
    for (i=0; i<4; i++) {
        file_phases[i] = (char*)malloc(100*sizeof(char));
    }
    
    //reading input file
    printf("Reading input file...\n");
    input(argc[1],&s_step,sub_const,&n0_sub,&n1_sub,&n2_sub,file_phases);
    //printf("%.3f %s %d %d %d\n",s_step,sub_const,n0_sub,n1_sub,n2_sub);
    //printf("%s\n",filename);
    printf("Done.\n\n");
    
    subtraction_constant(sub_const,&n0,&n1,&n2);
    //printf("%d %d %d\n",n0,n1,n2);
    printf("Dertermine the system for subtraction constant '%s'\n\n",sub_const);
    
    //import scattering phase shifts
    printf("Import scattering phase shifts...\n");
    sprintf(filename,"Input/Phases/%s",file_phases[0]);
    phase_input(filename,&delta0_spline,&L2_delta0,&delta0_const);
    sprintf(filename,"Input/Phases/%s",file_phases[1]);
    phase_input(filename,&delta1_spline,&L2_delta1,&delta1_const);
    sprintf(filename,"Input/Phases/%s",file_phases[2]);
    phase_input(filename,&delta2_spline,&L2_delta2,&delta2_const);
    sprintf(filename,"Input/Phases/%s",file_phases[3]);
    phase_input(filename,&delta_etapi_spline,&L2_delta_etapi,&delta_etapi_const);
    printf("Done.\n\n");
    
    //allocate gsl spline accelerators for the phase shifts
    acc_delta0 = gsl_interp_accel_alloc();
    acc_delta1 = gsl_interp_accel_alloc();
    acc_delta2 = gsl_interp_accel_alloc();
    acc_delta_etapi = gsl_interp_accel_alloc();
    
    //exit(0);
    
    printf("s_th = %.6f MPION**2\n",s0);
    printf("L2 = %.6f MPION**2\n\n",L2);
    
    s0 = s_pipi;
    L2 = 1000.;
    
    s_max = L2;
    s_min = integration_s_plus(s_max);
    
    eps = 1.0e-06;
    eps_b = 0.1; //treatment of dividing by 0 at b=METAP_P_MPION_SQUARED
    
    //for singularity in M_hat at a=(METAP-MPION)^2
    N_sin = (int)floor(log2(s_step/eps));
    //N_sin = 0;
    printf("N_singularity=%d %.2e\n\n",N_sin,s_step*pow(0.5,N_sin));
    
    s_step_low_high = 1.; //s_step for low and high energy ends
    //calculating number of grid points
    grid_point_size(N,N_sin,s_step,s_step_low_high,s_min,s_max,s0,eps_b);
    
    //number of grid points for angular averages
    N_M_avg[0] = N[2]; //s_th to 0.5*(METAP^2-MPION^2)
    N_M_avg[1] = N[3]+N[9]; //0.5*(METAP^2-MPION^2) to (META+MPION)^2 + (META+MPION)^2 to (METAP-MPION)^2
    N_M_avg[2] = N[4]; //(METAP-MPION)^2 to (METAP+MPION)^2
    N_M_avg[3] = N[5]+N[6]; //(METAP+MPION)^2 to 100*MPION^2 + 100*MPION^2 to s_max
    
    //number of grid points for Omnes-functions, amplitudes and inhomogeneities
    N_M_int[0] = N_M_avg[0]+N_M_avg[1]+N_M_avg[2]+N_M_avg[3]; //upper rim of the cut [s_th,L2]
    N_M_int[1] = N[7]; // lower rim of the cut [s_th,MPION*(METAP+MPION)]
    N_M_int[2] = N[0]+N[1]; //left from the cut [s_min,0]+[0,s_th]
    N_M_int[3] = N[8]; //complex plane
    N_M_int[4] = N[8]; //complex plane
    N_M_int[5] = N[8]; //complex plane
    N_M_int[6] = N[8]; //complex plane
    
    //number of grid points for tilde and hat functions
    N_M_hat[0] = N_M_avg[0]+N_M_avg[1]; //[s_th,a]
    N_M_hat[1] = N_M_avg[2]+N_M_avg[3]; //[a,L2]
    
    printf("Number of grid points for the Omnes functions, amplitudes and inhomogeneities:\n");
    for (i=0; i<4; i++) {
        printf("N_M_int(%d)=%d\n",i+1,N_M_int[i]);
    }
    printf("\n");
    
    printf("Number of grid points for angular averages:\n");
    for (i=0; i<4; i++) {
        printf("N_M_avg(%d)=%d\n",i+1,N_M_avg[i]);
    }
    printf("\n");
    
    printf("Number of grid points for tilde and hat functions:\n");
    for (i=0; i<2; i++) {
        printf("N_M_hat(%d)=%d\n",i+1,N_M_hat[i]);
    }
    printf("\n");
    
    s_below_cut = (double*)malloc(N_M_int[2]*sizeof(double));
    s_cv_plus = (double*)malloc(N_M_int[0]*sizeof(double));
    s_cv_minus = (double*)malloc(N_M_int[1]*sizeof(double));
    
    ya = (double*)malloc(N_M_int[3]*sizeof(double));
    yb = (double*)malloc(N_M_int[3]*sizeof(double));
    
    ys = (double**)malloc(2*sizeof(double*));
    for (i=0; i<2; i++) {
        ys[i] = (double*)malloc(N_M_int[3]*sizeof(double));
    }
    s_cmp = (complex**)malloc(2*sizeof(complex*));
    for (i=0; i<2; i++) {
        s_cmp[i] = (complex*)malloc(N_M_int[3]*sizeof(complex));
    }
    
    sa_complex = (complex*)malloc(N_M_int[3]*sizeof(complex));
    sb_complex = (complex*)malloc(N_M_int[3]*sizeof(complex));
    
    build_grid(s_cv_plus,s_cv_minus,s_below_cut,s_cmp,ys,s_min,s_max,s0,eps,eps_b,a,b,c,d,'0','+',N,N_sin);
    
    //allocate memory for omnes_I functions
    complex **omnes0 = (complex **)malloc(5*sizeof(complex *));
    for (i=0; i<5; i++) {
        omnes0[i] = (complex *)malloc(N_M_int[i]*sizeof(complex));
    }
    
    complex **omnes1 = (complex **)malloc(5*sizeof(complex *));
    for (i=0; i<5; i++) {
        omnes1[i] = (complex *)malloc(N_M_int[i]*sizeof(complex));
    }
    
    complex **omnes2 = (complex **)malloc(5*sizeof(complex *));
    for (i=0; i<5; i++) {
        omnes2[i] = (complex *)malloc(N_M_int[i]*sizeof(complex));
    }
    
    //allocate memory for the inhomogeneity integrals
    complex **M0_inhom = (complex **)malloc(7*sizeof(complex *));
    for (i=0; i<7; i++) {
        M0_inhom[i] = (complex *)malloc(N_M_int[i]*sizeof(complex));
    }
    
    complex **M1_inhom = (complex **)malloc(7*sizeof(complex *));
    for (i=0; i<7; i++) {
        M1_inhom[i] = (complex *)malloc(N_M_int[i]*sizeof(complex));
    }
    
    complex **M2_inhom = (complex **)malloc(7*sizeof(complex *));
    for (i=0; i<7; i++) {
        M2_inhom[i] = (complex *)malloc(N_M_int[i]*sizeof(complex));
    }
    
    //allocate memory for the angular averages
    n=2;
    complex **M0_avg = (complex **)malloc(n*sizeof(complex *));
    for (i=0; i<n; i++) {
        M0_avg[i] = (complex *)malloc(N_M_int[0]*sizeof(complex));
    }
    
    n=3;
    complex **M1_avg = (complex **)malloc(n*sizeof(complex *));
    for (i=0; i<n; i++) {
        M1_avg[i] = (complex *)malloc(N_M_int[0]*sizeof(complex));
    }
    
    n=2;
    complex **M2_avg = (complex **)malloc(n*sizeof(complex *));
    for (i=0; i<n; i++) {
        M2_avg[i] = (complex *)malloc(N_M_int[0]*sizeof(complex));
    }
    
    //allocate memory for the hat functions
    complex_spline *M0_hat = (complex_spline *)malloc(2*sizeof(complex_spline));
    complex_spline *M1_hat = (complex_spline *)malloc(2*sizeof(complex_spline));
    complex_spline *M2_hat = (complex_spline *)malloc(2*sizeof(complex_spline));
    
    complex_spline_alloc(M0_hat,2,N_M_hat);
    complex_spline_alloc(M1_hat,2,N_M_hat);
    complex_spline_alloc(M2_hat,2,N_M_hat);
    
    //allocate memory for the amplitudes
    complex_spline *M0 = (complex_spline *)malloc(7*sizeof(complex_spline));
    complex_spline *M1 = (complex_spline *)malloc(7*sizeof(complex_spline));
    complex_spline *M2 = (complex_spline *)malloc(7*sizeof(complex_spline));
    
    complex_spline_alloc(M0,7,N_M_int);
    complex_spline_alloc(M1,7,N_M_int);
    complex_spline_alloc(M2,7,N_M_int);
    
    printf("Computing Omnes functions...\n");
    printf("Isospin 0...\n");
    build_omnes(omnes0,delta0,s_pipi,L2,delta0_const,s_cv_plus,s_cv_minus,s_below_cut,s_cmp,N_M_int);
    printf("Done.\n");
    printf("Isospin 1...\n");
    build_omnes(omnes1,delta1,s_pipi,L2,delta1_const,s_cv_plus,s_cv_minus,s_below_cut,s_cmp,N_M_int);
    printf("Done.\n");
    printf("Isospin 2...\n");
    build_omnes(omnes2,delta2,s_pipi,L2,delta2_const,s_cv_plus,s_cv_minus,s_below_cut,s_cmp,N_M_int);
    printf("Done.\n\n");
    
    printf("Setting initial inhomogeneities to zero...\n");
    build_inhomogenity_init(M0_inhom,N_M_int);
    build_inhomogenity_init(M1_inhom,N_M_int);
    build_inhomogenity_init(M2_inhom,N_M_int);
    printf("Done.\n\n");
    
    printf("Computing Amplitudes...\n");
    build_amplitude(M0,omnes0,M0_inhom,s_cv_plus,s_cv_minus,s_below_cut,s_cmp,ys,s_pipi,L2,N_M_int,n0);
    build_amplitude(M1,omnes1,M1_inhom,s_cv_plus,s_cv_minus,s_below_cut,s_cmp,ys,s_pipi,L2,N_M_int,n1);
    build_amplitude(M2,omnes2,M2_inhom,s_cv_plus,s_cv_minus,s_below_cut,s_cmp,ys,s_pipi,L2,N_M_int,n2);
    printf("Done.\n\n");
    
    printf("Computing angular averages...\n");
    for (n=0; n<2; n++) {
        build_M_avg(M0,M0_avg[n],s_cv_plus,N_M_avg,a,b,c,d,s0,'0','+',n);
        build_M_avg(M1,M1_avg[n],s_cv_plus,N_M_avg,a,b,c,d,s0,'0','+',n);
        build_M_avg(M2,M2_avg[n],s_cv_plus,N_M_avg,a,b,c,d,s0,'0','+',n);
    }
    build_M_avg(M1,M1_avg[n],s_cv_plus,N_M_avg,a,b,c,d,s0,'0','+',n);
    printf("Done.\n\n");
    
    s = 50.0;
    sprintf(filename,"convergence_%s",sub_const);
    fout = fopen(filename,"w");
    gsl_interp_accel *acc = gsl_interp_accel_alloc();
    fprintf(fout,"%d %+.5e %+.5e %+.5e %+.5e %+.5e %+.5e\n",0,gsl_spline_eval(M0[0].re,s,acc),gsl_spline_eval(M0[0].im,s,acc),gsl_spline_eval(M1[0].re,s,acc),gsl_spline_eval(M1[0].im,s,acc),gsl_spline_eval(M2[0].re,s,acc),gsl_spline_eval(M2[0].im,s,acc));
    
    for (i=1; i<20; i++) {
        printf("Iteration %d:\n",i);
        
        printf("Building inhomogenity integrands...\n");
        build_inhomogenity_integrand(M0_hat,M1_hat,M2_hat,F0,G0,F1,G1,F2,G2,M0_avg,M1_avg,M2_avg,omnes0[0],omnes1[0],omnes2[0],s_cv_plus,s0,L2,N_M_avg,N_sin,n0_sub,n1_sub,n2_sub);
        printf("Done.\n\n");
        
        printf("Computing inhomogenities...\n");
        build_inhomogenity_sqrt(M0_inhom,M0_hat,F0,G0,s0,a,L2,s_cv_plus,s_cv_minus,s_below_cut,s_cmp,N_M_int,n0_sub);
        build_inhomogenity_sqrt_cubed(M1_inhom,M1_hat,F1,G1,s0,a,L2,s_cv_plus,s_cv_minus,s_below_cut,s_cmp,N_M_int,n1_sub);
        build_inhomogenity_sqrt(M2_inhom,M2_hat,F2,G2,s0,a,L2,s_cv_plus,s_cv_minus,s_below_cut,s_cmp,N_M_int,n2_sub);
        printf("Done.\n\n");
        
        printf("Computing Amplitudes...\n");
        build_amplitude(M0,omnes0,M0_inhom,s_cv_plus,s_cv_minus,s_below_cut,s_cmp,ys,s0,L2,N_M_int,n0);
        build_amplitude(M1,omnes1,M1_inhom,s_cv_plus,s_cv_minus,s_below_cut,s_cmp,ys,s0,L2,N_M_int,n1);
        build_amplitude(M2,omnes2,M2_inhom,s_cv_plus,s_cv_minus,s_below_cut,s_cmp,ys,s0,L2,N_M_int,n2);
        printf("Done.\n\n");
        
        printf("Computing angular averages...\n");
        for (n=0; n<2; n++) {
            build_M_avg(M0,M0_avg[n],s_cv_plus,N_M_avg,a,b,c,d,s0,'0','+',n);
            build_M_avg(M1,M1_avg[n],s_cv_plus,N_M_avg,a,b,c,d,s0,'0','+',n);
            build_M_avg(M2,M2_avg[n],s_cv_plus,N_M_avg,a,b,c,d,s0,'0','+',n);
        }
        build_M_avg(M1,M1_avg[n],s_cv_plus,N_M_avg,a,b,c,d,s0,'0','+',n);
        printf("Done.\n\n");
        
        s = 50.0;
        fprintf(fout,"%d %+.5e %+.5e %+.5e %+.5e %+.5e %+.5e\n",i,gsl_spline_eval(M0[0].re,s,acc),gsl_spline_eval(M0[0].im,s,acc),gsl_spline_eval(M1[0].re,s,acc),gsl_spline_eval(M1[0].im,s,acc),gsl_spline_eval(M2[0].re,s,acc),gsl_spline_eval(M2[0].im,s,acc));
    }
    fclose(fout);
    gsl_interp_accel_free(acc);
    
    printf("Writing output...\n");
    
    Omnes_output("OmnesFunctions",omnes0,omnes1,omnes2,s_below_cut,s_cv_plus,N_M_int[2],N_M_int[0]);
    
    sprintf(filename,"M_tilde_%s",sub_const);
    Inhomogenities_output(filename,M0_avg,M1_avg,M2_avg,s_cv_plus,N_M_int[0]);
    
    sprintf(filename,"M_int_%s",sub_const);
    Dispersion_integral_output(filename,M0_inhom,M1_inhom,M2_inhom,s_below_cut,s_cv_plus,N_M_int[2],N_M_int[0]);
    
    sprintf(filename,"M_%s",sub_const);
    Amplitudes_output(filename,M0,M1,M2,s_below_cut,s_cv_plus,N_M_int[2],N_M_int[0]);
    
    printf("Done.\n\n");
    
    free(s_cv_plus);
    free(s_cv_minus);
    free(s_below_cut);
    
    for (i=0; i<4; i++) {
        free(omnes0[i]);
    }
    free(omnes0);
    
    for (i=0; i<4; i++) {
        free(omnes1[i]);
    }
    free(omnes1);
    
    for (i=0; i<4; i++) {
        free(omnes2[i]);
    }
    free(omnes2);
    
    for (i=0; i<4; i++) {
        free(M0_inhom[i]);
    }
    free(M0_inhom);
    
    for (i=0; i<4; i++) {
        free(M1_inhom[i]);
    }
    free(M1_inhom);
    
    for (i=0; i<4; i++) {
        free(M2_inhom[i]);
    }
    free(M2_inhom);
    
    n=2;
    for (i=0; i<n; i++) {
        free(M0_avg[i]);
    }
    free(M0_avg);
    
    n=3;
    for (i=0; i<n; i++) {
        free(M1_avg[i]);
    }
    free(M1_avg);
    
    n=2;
    for (i=0; i<n; i++) {
        free(M2_avg[i]);
    }
    free(M2_avg);
    
    complex_spline_free(M0_hat,2);
    complex_spline_free(M1_hat,2);
    complex_spline_free(M2_hat,2);
    free(M0_hat);
    free(M1_hat);
    free(M2_hat);
    
    complex_spline_free(M0,4);
    complex_spline_free(M1,4);
    complex_spline_free(M2,4);
    free(M0);
    free(M1);
    free(M2);
    
    gsl_spline_free(delta0_spline);
    gsl_spline_free(delta1_spline);
    gsl_spline_free(delta2_spline);
    gsl_spline_free(delta_etapi_spline);
    
    gsl_interp_accel_free(acc_delta0);
    gsl_interp_accel_free(acc_delta1);
    gsl_interp_accel_free(acc_delta2);
    gsl_interp_accel_free(acc_delta_etapi);
    
    t = (clock()-t)/CLOCKS_PER_SEC;
    
    printf("\ncalculation time = %.3fs\n\n", t);
    
    return 0;
}

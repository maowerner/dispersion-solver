#
#-Wall -pedantic

cd src

gcc-7 $1\
    -c -std=c++14 \
    Basic.c PhaseShifts.c Grid.c Omnes.c Singularity.c AngularAverages.c InhomogenitySqrt.c InhomogenitySqrtCubed.c Amplitude.c InputOutput.c main_etap_eta2pi.c

gcc-7 -lm -O3 $(gsl-config --libs) *.o -o EtaPrime

rm *.o

mv EtaPrime ../EtaPrime

cd ..

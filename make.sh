#
#-Wall -pedantic

cd src

gcc $1\
    -c -std=c99 \
    Basic.c Grid.c Omnes.c Singularity.c AngularAverages.c InhomogenitySqrt.c InhomogenitySqrtCubed.c Amplitude.c InputOutput.c main.c

gcc -lm -O3 $(gsl-config --libs) *.o -o EtaPrime

rm *.o

mv EtaPrime ../EtaPrime

cd ..

#ifndef _PHASESHIFTS_H_
#define _PHASESHIFTS_H_

delta_params delta_high_energy_tale(delta_params delta);

//pipi I=0 S-wave scattering phase shift
double delta0(double s);

//pipi I=1 P-wave scattering phase shift
double delta1(double s);

//pipi I=2 S-wave scattering phase shift
double delta2(double s);

//etapi I=1 S-wave scattering phase shift
double delta_etapi(double s);

#endif

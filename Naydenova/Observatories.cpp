#include "Observatories.h"

/*
Code  Long.    cos       sin     Name
007   2.33675 0.659470 +0.749223 Paris
J80 359.1083  0.70862  +0.70323  Sainte Helene
080  28.9667  0.75566  +0.65278  Istanbul
724 260.8053  0.94388  +0.33026  National Observatory, Tacubaya
C13   9.10031 0.698332 +0.713430 Como
E12 149.0642  0.85563  -0.51621  Siding Spring Survey
327 117.5750  0.76278  +0.64470  Peking Observatory, Xinglong Station
Z19 342.11094 0.877701 +0.478380 La Palma-TNG
*/

Observatories::Observatories(double longitude, double cos, double sin): longitude(longitude), cos_l(cos), sin_l(sin) {
    vec = new double[3];
    CylindrCoordToCartesian();
};

void Observatories::CylindrCoordToCartesian() {
    vec[0] = R_CONST * cos(this->longitude) * cos_l;
    vec[1] = R_CONST * sin(this->longitude) * cos_l;
    vec[2] = R_CONST * sin_l;
}


double* Observatories::get_coords(){
    return vec;
}


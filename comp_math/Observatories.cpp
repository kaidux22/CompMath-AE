#include "Observatories.h"

/*
Code  Long.    cos       sin     Name
007   2.33675 0.659470 +0.749223 Paris
J80 359.1083  0.70862  +0.70323  Sainte Helene
080  28.9667  0.75566  +0.65278  Istanbul
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


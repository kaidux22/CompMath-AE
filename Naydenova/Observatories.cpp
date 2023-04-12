#include "Observatories.h"



Observatories::Observatories(double longitude, double cos, double sin): longitude(longitude), cos_l(cos), sin_l(sin) {
    vec = new double[3];
    this->angle = this->longitude * PI /180.0;
    SphericalCoordToCartesian();
};

void Observatories::SphericalCoordToCartesian() {
    vec[0] = R_CONST * cos(this->angle) * cos_l;
    vec[1] = R_CONST * sin(this->angle) * cos_l;
    vec[2] = R_CONST * sin_l;
}


double* Observatories::get_coords(){
    return vec;
}


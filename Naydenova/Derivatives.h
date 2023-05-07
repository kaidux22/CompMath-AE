#pragma once

#include <cmath>
#define MU_CONST 398600.4415 // км^3/с^2
#define J2 1.75553e10

double dx(double* vec) ;

double dy(double* vec) ;

double dz(double* vec) ;

double dg_dx(double* vec, double* coords);
double dg_dy(double* vec, double* coords);
double dg_dz(double* vec, double* coords);
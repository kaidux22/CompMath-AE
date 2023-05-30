#pragma once
#include "sofa/sofa.h"
#include "Converter.h"
#include <cmath>
#define MU_CONST 398600.4415 // км^3/с^2
#define J2 1.75553e10
using namespace std;

double dx(double* vec, double mu, double J) ;

double dy(double* vec, double mu, double J) ;

double dz(double* vec, double mu, double J) ;

double dg_dx(double* vec, double* coords);
double dg_dy(double* vec, double* coords);
double dg_dz(double* vec, double* coords);
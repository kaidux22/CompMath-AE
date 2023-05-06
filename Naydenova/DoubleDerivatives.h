#pragma once

#include <cmath>

double dux_dx(double* vec, double mu, double J);
double dux_dy(double* vec, double mu, double J);
double dux_dz(double* vec, double mu, double J);
double duy_dx(double* vec, double mu, double J);
double duy_dy(double* vec, double mu, double J);
double duy_dz(double* vec, double mu, double J);
double duz_dx(double* vec, double mu, double J);
double duz_dy(double* vec, double mu, double J);
double duz_dz(double* vec, double mu, double J);
double dux_dmu(double* vec);
double duy_dmu(double* vec);
double duz_dmu(double* vec);
double dux_dJ(double* vec);
double duy_dJ(double* vec);
double duz_dJ(double* vec);
#include <iostream>
#include "GravPot.h"
#include "AdvancedFunction.h"
#define GENERAL_TIME 86400

using namespace std;

void DormandPrince(double UTC, double h, const int N, double J, double mu, double* vec, double a[7][7], double b[7], double** k, double c[7], void (*f)(double*, double*, double, double, double));

double** integrate(double JD, double h, const int N, double* vec);

double** integrate_for_inverse(double JD, double h, const int N, double* vec, double J, double mu);
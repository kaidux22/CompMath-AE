#include <iostream>
#include "GravPot.h"
#define GENERAL_TIME 86400

using namespace std;

void DormandPrince(double UTC, double h, const int N, double* vec, double a[7][7], double b[7], double** k, void (*f)(double*, double));

double** intergrate(double UTC_start, double h, const int N, double* vec);
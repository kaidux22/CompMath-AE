#pragma once
#include <cmath>
#include <cassert>
#include <iostream>
#include "Vnm.h"

#define N_CONST 2
#define NU_CONST 398600.4415 // км^3/с^2
#define R_CONST 6378.1363 // км
#define MAX_ORD 2

using namespace std;

double GravPot(double *vec, ComplexNum(*d)(LegFunc&, int, int, double*)) {
	assert(vec[1] != 0);
	
	int N = N_CONST;
	double nu = NU_CONST, R = R_CONST;
	
	double J0 = 1.0; double J2 = -0.1082635854e-2;
	
	LegFunc Pmn = LegFunc(N + MAX_ORD, N + MAX_ORD, sqrt(vec[0] * vec[0] + vec[1] * vec[1]) / sqrt(vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2]));
	ComplexNum res = ComplexNum(0, 0);

	res = (ComplexNum)J0 * d(Pmn, 0, 0, vec) + (ComplexNum)pow(R, 2) * (ComplexNum)J2 * d(Pmn, 2, 0, vec);
	
	return res.Real() * (-nu);

}

double* GradV(double* vec, double time) {
	double* grad = new double[3];
	grad[0] = GravPot(vec, Vdx);
	grad[1] = GravPot(vec, Vdy);
	grad[2] = GravPot(vec, Vdz);

	return grad;
}

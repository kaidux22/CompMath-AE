#pragma once
#include <cmath>
#include <cassert>
#include <iostream>
#include "Vnm.h"

#define N_CONST 2
#define NU_CONST 398600.4415 // κμ^3/ρ^2
#define R_CONST 6378.1363 // κμ
#define MAX_ORD 2

using namespace std;

void GravPot(double *vec, ComplexNum(*dx)(LegFunc&, int, int, double*), ComplexNum(*dy)(LegFunc&, int, int, double* ), ComplexNum(*dz)(LegFunc&, int, int, double*)) {
	assert(vec[1] != 0);
	
	int N = N_CONST;
	double nu = NU_CONST, R = R_CONST;
	
	double J0 = 1.0; double J2 = -0.1082635854e-2;
	
	LegFunc Pmn = LegFunc(N + MAX_ORD, N + MAX_ORD, sqrt(vec[0] * vec[0] + vec[1] * vec[1]) / sqrt(vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2]));
	ComplexNum res1 = ComplexNum(0, 0);
	ComplexNum res2 = ComplexNum(0, 0);
	ComplexNum res3 = ComplexNum(0, 0);

	res1 = (ComplexNum)J0 * dx(Pmn, 0, 0, vec) + (ComplexNum)pow(R, 2) * (ComplexNum)J2 * dx(Pmn, 2, 0, vec);
	res2 = (ComplexNum)J0 * dy(Pmn, 0, 0, vec) + (ComplexNum)pow(R, 2) * (ComplexNum)J2 * dy(Pmn, 2, 0, vec);
	res3 = (ComplexNum)J0 * dz(Pmn, 0, 0, vec) + (ComplexNum)pow(R, 2) * (ComplexNum)J2 * dz(Pmn, 2, 0, vec);

	vec[0] = res1.Real() * (-nu);
	vec[1] = res2.Real() * (-nu);
	vec[2] = res3.Real() * (-nu);

}
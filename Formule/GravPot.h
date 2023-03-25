#ifndef COMPMATH_GRAVPOT
#define COMPMATH_GRAVPOT

#include <cmath>
#include <cassert>
#include <iostream>

#include "LegFunc.h"
#include "ComplexNums.h"

#define N_CONST 4
#define NU_CONST 398600.4415 // км^3/с^2
#define R_CONST 6378.1363 // км

using namespace std;

//функция для дальнейшего подсчёта градиента реккурентно
ComplexNum Vnm(LegFunc& Pmn, int n, int m, int x, int y, int z) {
	return (ComplexNum)(Pmn.ExtractValue(n, m) / (pow(sqrt(x * x + y * y + z * z), (double)(n + m + 1)))) * ComplexNum(x, y).Pow(m);
}

double GravPot(double x, double y, double z, ComplexNum(*func)(LegFunc&, int, int, int, int, int)) {
	assert(y != 0);
	int N = N_CONST;
	double nu = NU_CONST, R = R_CONST;
	//первый индекс - m, второй индекс - n
	//С[0][n] = Jn
	double Cmn[5][5] = { { 1.0, 0.0, -0.1082635854e-2, 0.2532435346e-5, 0.1619331205e-5 },
					   {0.0, 0.0, -0.3504890360e-9, 0.2192798802e-5, -0.5087253036e-6},
					   {0.0, 0.0, 0.1574536043e-5, 0.3090160446e-6, 0.7841223074e-7},
					   {0.0, 0.0, 0.0, 0.1005588574e-6, 0.5921574319e-7},
					   {0.0, 0.0, 0.0, 0.0, -0.3982395740e-8} };
	// первый индекс - m, второй индекс - n
	double Smn[5][5] = { {0.0, 0.0, 0.0, 0.0, 0.0},
					   {0.0, 0.0, 0.1635406077e-8, 0.2680118938e-6, -0.4494599352e-6},
					   {0.0, 0.0, -0.9038680729e-6, -0.2114023978e-6, 0.1481554569e-6},
					   {0.0, 0.0, 0.0, 0.1972013239e-6, -0.1201129183e-7},
					   {0.0, 0.0, 0.0, 0.0, 0.6525605810e-8} };

	LegFunc Pmn = LegFunc(N, N, sqrt(x * x + y * y) / sqrt(x * x + y * y + z * z));
	Pmn.PrintMaxtrix();
	ComplexNum res = ComplexNum(0, 0);

	for (int n = 0; n <= N; n++) {
		if (n == 1) {
			continue;
		}
		for (int m = 0; m <= n; m++) {
			res = res + (ComplexNum)pow(R, n) * ComplexNum(Cmn[m][n], -Smn[m][n]) * func(Pmn, n, m, x, y, z);
		}
	}

	return res.Real() * (-nu);
}

#endif //COMPMATH_GRAVPOT

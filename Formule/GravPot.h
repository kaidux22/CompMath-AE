#ifndef COMPMATH_GRAVPOT
#define COMPMATH_GRAVPOT

#include <cmath>
#include <cassert>
#include <iostream>

#include "LegFunc.h"
#include "ComplexNums.h"
#include "Vnm.h"
#include "Converter.h"

#define N_CONST 2
#define NU_CONST 398600.4415 // км^3/с^2
#define R_CONST 6378.1363 // км
#define MAX_ORD 2 //наибольшая степень производной

using namespace std;

double GravPot(double* vec, ComplexNum(*func)(LegFunc&, int, int, double*)) {
	assert(vec[1] != 0);
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

	//Создаю таблицу значений полиномов Лежандра до N + 2 степени и порядка
	LegFunc Pmn = LegFunc(N + MAX_ORD, N + MAX_ORD, sqrt(vec[0] * vec[0] + vec[1] * vec[1]) / sqrt(vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2]));
	//Pmn.PrintMaxtrix();
	ComplexNum res = ComplexNum(0, 0);


	for (int n = 0; n <= N; n++) {
		if (n == 1) {
			continue;
		}
		for (int m = 0; m <= n; m++) {
			res = res + (ComplexNum)pow(R, n) * ComplexNum(Cmn[m][n], -Smn[m][n]) * func(Pmn, n, m, vec);
		}
	}

	return res.Real() * (-nu);
}

//метод возвращается градиент гравитационного потенциала
double* GradV(double* vec, double UTC) {
	double rotateMatrix[3][3];

	iauC2t06a(UTC + 37 + 32.184, 0, UTC, 0, 0, 0, rotateMatrix);
	
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			cout << rotateMatrix[i][j] << ' ';
		}
		cout << endl;
	}

	BodySpaceFixed(vec, UTC, rotateMatrix);

	double* grad = new double[3];
	grad[0] = GravPot(vec, Vdx);
	grad[1] = GravPot(vec, Vdy);
	grad[2] = GravPot(vec, Vdz);

	Trans(rotateMatrix);

	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			cout << rotateMatrix[i][j] << ' ';
		}
		cout << endl;
	}

	BodySpaceFixed(vec, UTC, rotateMatrix);

	return grad;
}

#endif //COMPMATH_GRAVPOT
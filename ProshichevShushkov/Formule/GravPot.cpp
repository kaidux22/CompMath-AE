#include "GravPot.h"

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

/*
y'(t) = v(t)
v'(t) = f
 */

void DormandPrince(double t, double h, const int N, double* vec, double a[7][7], double b[7], int integrate_numder, double* (*f)(double* vec, double)) {

	double** k = new double* [7];
	for (int i = 0; i < 7; i++) {
		k[i] = new double[N];
		for (int j = 0; j < N; j++) {
			k[i][j] = 0;
		}
	}

	for (int i = 0; i < 7; i++) {
		for (int j = 0; j < N; j++) {
			k[i][j] = vec[j];
			for (int t = 0; t < 7; t++) {
				k[i][j] += a[i][t] * h * k[t][j];
			}
		}
		if (!integrate_numder) {
			k[i] = f(k[i], t);
		}

	}

	for (int i = 0; i < N; i++) {
		for (int j = 0; j < 7; j++) {
			vec[i] += b[j] * k[j][i] * h;
		}

	}

}

void intergrate(double UTC_start, double h, const int N, double* vec) {
	
	double a[7][7] = { {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
					   {1.0 / 5.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
					   {3.0 / 40.0, 9.0 / 40.0, 0.0, 0.0, 0.0, 0.0, 0.0},
					   {44.0 / 45.0, -56.0 / 15.0, 32.0 / 9.0, 0.0, 0.0, 0.0, 0.0 },
					   {19372.0 / 6561.0, -25360.0 / 2187.0, 64448.0 / 6561.0, -212.0 / 729.0, 0.0, 0.0, 0.0 },
					   {9017.0 / 3168.0, -355.0 / 33.0, 46732.0 / 5247.0, 49.0 / 176.0, -5103.0 / 18656.0, 0.0, 0.0},
					   {35.0 / 384.0, 0.0, 500.0 / 1113.0, 125.0 / 192.0, -2187.0 / 6784.0, 11.0 / 84.0, 0.0},
	};

	double b[7] = { 5179.0 / 57600.0, 0.0, 7571.0 / 16695.0, 393.0 / 640.0, -92097.0 / 339200.0, 187.0 / 2100.0, 1.0 / 40.0 };
	
	double UTC = 2451545.0 * 24.0 * 60.0 * 60.0;
	
	// 86400 секунд в сутках
	for (int i = 0; i < 86400 / h; i++) {
		for (int i = 0; i < 2; i++) {
			DormandPrince(UTC, h, N, vec, a, b, i, GradV);
			
		}
		UTC += h;
	}
}
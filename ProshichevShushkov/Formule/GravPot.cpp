#include "GravPot.h"



double GravPot(double* vec, ComplexNum(*func)(LegFunc&, int, int, double*)) {
	int N = N_CONST;
	double R = R_CONST;
	//первый индекс - m, второй индекс - n
	//С[0][n] = Jn
	double Cmn[5][5] = { { -1.0, 0.0, 0.1082635854e-2, -0.2532435346e-5, -0.1619331205e-5 },
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
	LegFunc Pmn = LegFunc(N + MAX_ORD, N + MAX_ORD, vec[2] / sqrt(vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2]));
	//Pmn.PrintMaxtrix();
	ComplexNum res = ComplexNum(0, 0);


	for (int n = 0; n <= N; n++) {
		if (n == 1) {
			continue;
		}
		for (int m = 0; m <= n; m++) {
			res = res + (ComplexNum(R, 0).Pow(n) * ComplexNum(Cmn[m][n], -Smn[m][n])) * func(Pmn, n, m, vec);
		}
	}

	return res.Real();
}

double DerivativedVdC(double* vec, int n, int m, ComplexNum(*func)(LegFunc&, int, int, double*)){
    int N = N_CONST;
    LegFunc Pmn = LegFunc(N + MAX_ORD, N + MAX_ORD, vec[2] / sqrt(vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2]));
    ComplexNum res = ComplexNum(R_CONST,0).Pow(n) * func(Pmn, n, m, vec);
    return res.Real();
}

double DerivativedVdS(double* vec, int n, int m, ComplexNum(*func)(LegFunc&, int, int, double*)){
    int N = N_CONST;
    LegFunc Pmn = LegFunc(N + MAX_ORD, N + MAX_ORD, vec[2] / sqrt(vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2]));
    ComplexNum res = ComplexNum(0,-1) * ComplexNum(R_CONST,0).Pow(n) * func(Pmn, n, m, vec);
    return res.Real();
}

double DerivativedVdM(double* vec, ComplexNum(*func)(LegFunc&, int, int, double*)){
    return (NU_CONST/M_CONST) * GravPot(vec, func);
}

//метод возвращается градиент гравитационного потенциала
void GradV(double* x, double* vec, double JD) {

    vec[0] = x[3];
    vec[1] = x[4];
    vec[2] = x[5];
    vec[3] = 0;
    vec[4] = 0;
    vec[5] = 0;

    double rotateMatrix[3][3];

    iauC2t06a(JD + (37.0 + 32.184) / 86400.0, 0, JD, 0, 0, 0, rotateMatrix);

    changeCoords(rotateMatrix, x, 0);

    double *grad = new double[3];
    grad[0] = -NU_CONST * GravPot(x, Vdx);
    grad[1] = -NU_CONST * GravPot(x, Vdy);
    grad[2] = -NU_CONST * GravPot(x, Vdz);

    Transposition(rotateMatrix);

    changeCoords(rotateMatrix, grad, 0);

    for (int i = 0; i < 3; i++) {
        vec[i + 3] = grad[i];
    }
}

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



void DormandPrince(double JD, double h, const int N, double* vec, double a[7][7], double b[7], double** k, double c[7], void (*f)(double*, double*, double)) {

    double* x = new double[N];

    for (int i = 0 ; i < 7; i++) {
        for (int j = 0; j < N; j++) {
            x[j] = vec[j];
            for (int t = 0; t < i; t++) {
                x[j] += a[i][t] * h * k[t][j];
            }
        }
        f(x, k[i], JD + c[i]* h);
    }

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < 7; j++) {
            vec[i] += b[j] * k[j][i] *h ;
        }
    }

    delete[] x;

}

double** intergrate(double JD, double h, const int N, double* vec) {

    double a[7][7] = { {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
                       {1.0 / 5.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
                       {3.0 / 40.0, 9.0 / 40.0, 0.0, 0.0, 0.0, 0.0, 0.0},
                       {44.0 / 45.0, -56.0 / 15.0, 32.0 / 9.0, 0.0, 0.0, 0.0, 0.0 },
                       {19372.0 / 6561.0, -25360.0 / 2187.0, 64448.0 / 6561.0, -212.0 / 729.0, 0.0, 0.0, 0.0 },
                       {9017.0 / 3168.0, -355.0 / 33.0, 46732.0 / 5247.0, 49.0 / 176.0, -5103.0 / 18656.0, 0.0, 0.0},
                       {35.0 / 384.0, 0.0, 500.0 / 1113.0, 125.0 / 192.0, -2187.0 / 6784.0, 11.0 / 84.0, 0.0},
    };
    double c[7] = {0.0, 1.0 / 5.0, 3.0/10.0, 4.0/5.0, 8.0/9.0, 1.0, 1.0 };

    double b[7] = { 5179.0 / 57600.0, 0.0, 7571.0 / 16695.0, 393.0 / 640.0, -92097.0 / 339200.0, 187.0 / 2100.0, 1.0 / 40.0 };

    double** k = new double* [7];
    for (int i = 0; i < 7; i++) {
        k[i] = new double[N];
        for (int j = 0; j < N; j++) {
            k[i][j] = 0;
        }
    }

    int cnt = GENERAL_TIME / h;
    double** orbit = new double* [cnt];
    // 86400 секунд в сутках
    for (int i = 0; i < cnt; i++) {
        DormandPrince(JD, h, N, vec, a, b, k, c, GradV);
        orbit[i] = new double[7];
        orbit[i][0] = JD, orbit[i][1] = vec[0], orbit[i][2] = vec[1], orbit[i][3] = vec[2], orbit[i][4] = vec[3], orbit[i][5] = vec[4], orbit[i][6] = vec[5];
        JD += h / 86400.0;
    }

    for (int i = 0; i < N; i++) {
        delete[] k[i];
    }
    delete[] k;

    return orbit;

}

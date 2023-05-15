#include "RightPart.h"

double GravPotWithParams(double* vec, Matrix<double> *params, ComplexNum(*func)(LegFunc&, int, int, double*)) {

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
	int cnt = 13;
    //Cmn
    for(int n = 2; n < 5; n++){
        for(int m = 0; m <= n; m++){
            Cmn[m][n] = params->Get(cnt, 0);
            cnt++;
        }
    }

    //Smn
    for(int n = 2; n < 5; n++){
        for(int m = 1; m <= n; m++){
            Smn[m][n] = params->Get(cnt, 0);
            cnt++;
        }
    }

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

	return params->Get(12, 0) * res.Real();
}

Matrix<double> *MatrixdFdX(double *x, Matrix<double> *params){
    Matrix<double> *dFdX = new Matrix<double>(12, 12);
    for(int i = 0; i < 6; i++){
        dFdX->Set(i, i + 6, 1);
    }

    dFdX->Set(6, 0, -GravPotWithParams(x, params, Vdxdx)), dFdX->Set(6, 1, -GravPotWithParams(x, params, Vdxdy)), dFdX->Set(6, 2, -GravPotWithParams(x, params, Vdxdz));
	dFdX->Set(6, 3, -GravPotWithParams(x, params, Vdxdx)), dFdX->Set(6, 4, -GravPotWithParams(x, params, Vdxdy)), dFdX->Set(6, 5, -GravPotWithParams(x, params, Vdxdz));

    dFdX->Set(7, 0, -GravPotWithParams(x, params, Vdxdy)), dFdX->Set(7, 1, -GravPotWithParams(x, params, Vdydy)), dFdX->Set(7, 2, -GravPotWithParams(x, params, Vdydz));
	dFdX->Set(7, 3, -GravPotWithParams(x, params, Vdxdy)), dFdX->Set(7, 4, -GravPotWithParams(x, params, Vdydy)), dFdX->Set(7, 5, -GravPotWithParams(x, params, Vdydz));

    dFdX->Set(8, 0, -GravPotWithParams(x, params, Vdxdz)), dFdX->Set(8, 1, -GravPotWithParams(x, params, Vdydz)), dFdX->Set(8, 2, -GravPotWithParams(x, params, Vdzdz));
	dFdX->Set(8, 3, -GravPotWithParams(x, params, Vdxdz)), dFdX->Set(8, 4, -GravPotWithParams(x, params, Vdydz)), dFdX->Set(8, 5, -GravPotWithParams(x, params, Vdzdz));

	dFdX->Set(9, 0, -GravPotWithParams(x + 3, params, Vdxdx)), dFdX->Set(9, 1, -GravPotWithParams(x + 3, params, Vdxdy)), dFdX->Set(9, 2, -GravPotWithParams(x + 3, params, Vdxdz));
	dFdX->Set(9, 3, -GravPotWithParams(x + 3, params, Vdxdx)), dFdX->Set(9, 4, -GravPotWithParams(x + 3, params, Vdxdy)), dFdX->Set(9, 5, -GravPotWithParams(x + 3, params, Vdxdz));

    dFdX->Set(10, 0, -GravPotWithParams(x + 3, params, Vdxdy)), dFdX->Set(10, 1, -GravPotWithParams(x + 3, params, Vdydy)), dFdX->Set(10, 2, -GravPotWithParams(x + 3, params, Vdydz));
	dFdX->Set(10, 3, -GravPotWithParams(x + 3, params, Vdxdy)), dFdX->Set(10, 4, -GravPotWithParams(x + 3, params, Vdydy)), dFdX->Set(10, 5, -GravPotWithParams(x + 3, params, Vdydz));

    dFdX->Set(11, 0, -GravPotWithParams(x + 3, params, Vdxdz)), dFdX->Set(11, 1, -GravPotWithParams(x + 3, params, Vdydz)), dFdX->Set(11, 2, -GravPotWithParams(x + 3, params, Vdzdz));
	dFdX->Set(11, 3, -GravPotWithParams(x + 3, params, Vdxdz)), dFdX->Set(11, 4, -GravPotWithParams(x + 3, params, Vdydz)), dFdX->Set(11, 5, -GravPotWithParams(x + 3, params, Vdzdz));

	return dFdX;
}

Matrix<double> *MatrixdFdParam(double *x, Matrix<double> *params){
	Matrix<double> *dFdParam = new Matrix<double>(12, 34);

    dFdParam->Set(6, 12, -DerivativedVdGM(x, params, Vdx)), dFdParam->Set(7, 12, -DerivativedVdGM(x, params, Vdy)), dFdParam->Set(8, 12, -DerivativedVdGM(x, params, Vdz));
	dFdParam->Set(9, 12, -DerivativedVdGM(x + 3, params, Vdx)), dFdParam->Set(10, 12, -DerivativedVdGM(x + 3, params, Vdy)), dFdParam->Set(11, 12, -DerivativedVdGM(x + 3, params, Vdz));

	int cnt = 13;
    //Cmn
    for(int n = 2; n < 5; n++){
        for(int m = 0; m <= n; m++){
            dFdParam->Set(6, cnt, -DerivativedVdC(x, params, n, m, Vdx)), dFdParam->Set(7, cnt, -DerivativedVdC(x, params, n, m, Vdy)), dFdParam->Set(8, cnt, -DerivativedVdC(x, params, n, m, Vdz));
			dFdParam->Set(9, cnt, -DerivativedVdC(x + 3, params, n, m, Vdx)), dFdParam->Set(10, cnt, -DerivativedVdC(x + 3, params, n, m, Vdy)), dFdParam->Set(11, cnt, -DerivativedVdC(x + 3, params, n, m, Vdz));
            cnt++;
        }
    }

    //Smn
    for(int n = 2; n < 5; n++){
        for(int m = 1; m <= n; m++){
            dFdParam->Set(6, cnt, -DerivativedVdS(x, params, n, m, Vdx)), dFdParam->Set(7, cnt, -DerivativedVdS(x, params, n, m, Vdy)), dFdParam->Set(8, cnt, -DerivativedVdS(x, params, n, m, Vdz));
			dFdParam->Set(9, cnt, -DerivativedVdS(x + 3, params, n, m, Vdx)), dFdParam->Set(10, cnt, -DerivativedVdS(x + 3, params, n, m, Vdy)), dFdParam->Set(11, cnt, -DerivativedVdS(x + 3, params, n, m, Vdz));
            cnt++;
        }
    }

    return dFdParam;
}

void RightPart(double* x, double* vec, double JD, Matrix<double> *params) {
    vec[0] = x[6];
    vec[1] = x[7];
    vec[2] = x[8];
    vec[3] = x[9];
    vec[4] = x[10];
    vec[5] = x[11];

	vec[6] = 0;
	vec[7] = 0;
	vec[8] = 0;
	vec[9] = 0;
	vec[10] = 0;
	vec[11] = 0;

	/*
	for(int i = 0; i < 12; i++){
		cout << vec[i] << endl;
	}
	cout << endl;
	*/

    double rotateMatrix[3][3];

    iauC2t06a(JD + (37.0 + 32.184) / 86400.0, 0, JD, 0, 0, 0, rotateMatrix);

    for(int i = 0; i <  420 / 3; i++)
        changeCoords(rotateMatrix, x, 3 * i);

    double *grad = new double[6];
	
    grad[0] = -GravPotWithParams(x, params, Vdx);
    grad[1] = -GravPotWithParams(x, params, Vdy);
    grad[2] = -GravPotWithParams(x, params, Vdz);
	
	grad[3] = -GravPotWithParams(x + 3, params, Vdx);
	grad[4] = -GravPotWithParams(x + 3, params, Vdy);
	grad[5] = -GravPotWithParams(x + 3, params, Vdz);
	
	Matrix<double> *dFdX = MatrixdFdX(x, params);
	Matrix<double> *dXdParam = new Matrix<double>(x + 12, 12, 34);
    Matrix<double> *dFdParam = MatrixdFdParam(x, params);

	double *res = (*dFdX * *dXdParam + *dFdParam).TransToVector();
	for(int i = 12; i < 12 * 34; i++){
		vec[i] = res[i - 12];
	}

	delete dFdX;
	delete dXdParam;
	delete dFdParam;

    Transposition(rotateMatrix);

    for (int i = 0; i < 6; i++) {
        vec[i + 6] = grad[i];
    }

	for(int i = 2; i < 12 * 35 / 3; i++){
		changeCoords(rotateMatrix, vec, 3 * i);
	}

}

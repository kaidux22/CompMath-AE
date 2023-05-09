#include "RightPart.h"

double GravPotWithParams(double* vec, Matrix<double> *params, ComplexNum(*func)(LegFunc&, int, int, double*)) {
	int N = N_CONST;
	double R = R_CONST;
	//первый индекс - m, второй индекс - n
	//С[0][n] = Jn
	double Cmn[5][5] = { { -1.0, 0.0, params->Get(13, 0), params->Get(14, 0), params->Get(15, 0) },
					   {0.0, 0.0, params->Get(16, 0), params->Get(20, 0), params->Get(26, 0)},
					   {0.0, 0.0, params->Get(18, 0), params->Get(22, 0), params->Get(28, 0)},
					   {0.0, 0.0, 0.0, params->Get(24, 0), params->Get(30, 0)},
					   {0.0, 0.0, 0.0, 0.0, params->Get(32, 0)} };
	// первый индекс - m, второй индекс - n
	double Smn[5][5] = { {0.0, 0.0, 0.0, 0.0, 0.0},
					   {0.0, 0.0, params->Get(17, 0), params->Get(21, 0), params->Get(27, 0)},
					   {0.0, 0.0, params->Get(19, 0), params->Get(23, 0), params->Get(29, 0)},
					   {0.0, 0.0, 0.0, params->Get(25, 0), params->Get(31, 0)},
					   {0.0, 0.0, 0.0, 0.0, params->Get(33, 0)} };

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
    dFdX->Set(7, 0, -GravPotWithParams(x, params, Vdxdy)), dFdX->Set(7, 1, -GravPotWithParams(x, params, Vdydy)), dFdX->Set(7, 2, -GravPotWithParams(x, params, Vdydz));
    dFdX->Set(8, 0, -GravPotWithParams(x, params, Vdxdz)), dFdX->Set(8, 1, -GravPotWithParams(x, params, Vdydz)), dFdX->Set(8, 2, -GravPotWithParams(x, params, Vdzdz));

	dFdX->Set(9, 0, -GravPotWithParams(x + 6, params, Vdxdx)), dFdX->Set(9, 1, -GravPotWithParams(x + 6, params, Vdxdy)), dFdX->Set(9, 2, -GravPotWithParams(x + 6, params, Vdxdz));
    dFdX->Set(10, 0, -GravPotWithParams(x + 6, params, Vdxdy)), dFdX->Set(10, 1, -GravPotWithParams(x + 6, params, Vdydy)), dFdX->Set(10, 2, -GravPotWithParams(x + 6, params, Vdydz));
    dFdX->Set(11, 0, -GravPotWithParams(x + 6, params, Vdxdz)), dFdX->Set(11, 1, -GravPotWithParams(x + 6, params, Vdydz)), dFdX->Set(11, 2, -GravPotWithParams(x + 6, params, Vdzdz));	

	return dFdX;
}

Matrix<double> *MatrixdFdParam(double *x, Matrix<double> *params){
	Matrix<double> *dFdParam = new Matrix<double>(12, 34);

    dFdParam->Set(6, 12, -DerivativedVdGM(x, params, Vdx)), dFdParam->Set(7, 12, -DerivativedVdGM(x, params, Vdy)), dFdParam->Set(8, 12, -DerivativedVdGM(x, params, Vdz));
	dFdParam->Set(9, 12, -DerivativedVdGM(x + 12, params, Vdx)), dFdParam->Set(10, 12, -DerivativedVdGM(x + 12, params, Vdy)), dFdParam->Set(11, 12, -DerivativedVdGM(x + 6, params, Vdz));

	for(int i = 0; i < 3; i++){
		dFdParam->Set(6, 13 + i, -DerivativedVdC(x, params, 2 + i, 0, Vdx)), dFdParam->Set(7, 13 + i, -DerivativedVdC(x, params, 2 + i, 0, Vdy)), dFdParam->Set(8, 13 + i, -DerivativedVdC(x, params, 2 + i, 0, Vdz));
		dFdParam->Set(9, 13 + i, -DerivativedVdC(x + 12, params, 2 + i, 0, Vdx)), dFdParam->Set(10, 13 + i, -DerivativedVdC(x + 12, params, 2 + i, 0, Vdy)), dFdParam->Set(11, 13 + i, -DerivativedVdC(x + 12, params, 2 + i, 0, Vdz));
	}

	
	int cnt = 16;
	
	for(int n = 2; n < 5; n++){
        for(int m = 1; m <= n; m++){
			dFdParam->Set(6, cnt, -DerivativedVdC(x, params, n, m, Vdx)), dFdParam->Set(7, cnt, -DerivativedVdC(x, params, n, m, Vdy)), dFdParam->Set(8, cnt, -DerivativedVdC(x, params, n, m, Vdz));
			dFdParam->Set(9, cnt, -DerivativedVdC(x + 12, params, n, m, Vdx)), dFdParam->Set(10, cnt, -DerivativedVdC(x + 12, params, n, m, Vdy)), dFdParam->Set(11, cnt, -DerivativedVdC(x + 12, params, n, m, Vdz));

			dFdParam->Set(6, cnt + 1, -DerivativedVdS(x, params, n, m, Vdx)), dFdParam->Set(7, cnt + 1, -DerivativedVdS(x, params, n, m, Vdy)), dFdParam->Set(8, cnt + 1, -DerivativedVdS(x + 12, params, n, m, Vdz));
			dFdParam->Set(9, cnt + 1, -DerivativedVdS(x + 12, params, n, m, Vdx)), dFdParam->Set(10, cnt + 1, -DerivativedVdS(x + 12, params, n, m, Vdy)), dFdParam->Set(11, cnt + 1, -DerivativedVdS(x + 12, params, n, m, Vdz));
            cnt += 2;
        }
    }

    return dFdParam;
}

void RightPart(double* x, double* vec, double JD, Matrix<double> *params) {
    vec[0] = x[3];
    vec[1] = x[4];
    vec[2] = x[5];
    vec[3] = 0;
    vec[4] = 0;
    vec[5] = 0;

	vec[6] = x[9];
	vec[7] = x[10];
	vec[8] = x[11];
	vec[9] = 0;
	vec[10] = 0;
	vec[11] = 0;

    double rotateMatrix[3][3];

    iauC2t06a(JD + (37.0 + 32.184) / 86400.0, 0, JD, 0, 0, 0, rotateMatrix);

    for(int i = 0; i <  420 / 3; i++)
        changeCoords(rotateMatrix, x, 3 * i);

    double *grad = new double[6];
    grad[0] = -GravPotWithParams(x, params, Vdx);
    grad[1] = -GravPotWithParams(x, params, Vdy);
    grad[2] = -GravPotWithParams(x, params, Vdz);
	
	grad[3] = -GravPotWithParams(x + 6, params, Vdx);
	grad[4] = -GravPotWithParams(x + 6, params, Vdy);
	grad[5] = -GravPotWithParams(x + 6, params, Vdz);

	Matrix<double> *dFdX = MatrixdFdX(x, params);
	Matrix<double> *dXdParam = new Matrix<double>(vec + 12, 12, 34);
    Matrix<double> *dFdParam = MatrixdFdParam(x, params);

	double *res = (*dFdX * *dXdParam + *dFdParam).TransToVector();
	for(int i = 12; i < 12 * 34; i++){
		vec[i] = res[i - 12];
	}

	delete dFdX;
	delete dXdParam;
	delete dFdParam;

    Transposition(rotateMatrix);

    for (int i = 0; i < 3; i++) {
        vec[i + 3] = grad[i];
    }

	for(int i = 0; i < 12 * 34 / 3; i++){
		changeCoords(rotateMatrix, vec, 3 * i);
	}
}

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
    Matrix<double> *dFdX = new Matrix<double>(6, 6);
    for(int i = 0; i < 3; i++){
        dFdX->Set(i, i + 3, 1);
    }

    dFdX->Set(3, 0, -GravPotWithParams(x, params, Vdxdx)), dFdX->Set(3, 1, -GravPotWithParams(x, params, Vdxdy)), dFdX->Set(3, 2, -GravPotWithParams(x, params, Vdxdz));
    dFdX->Set(4, 0, -GravPotWithParams(x, params, Vdxdy)), dFdX->Set(4, 1, -GravPotWithParams(x, params, Vdydy)), dFdX->Set(4, 2, -GravPotWithParams(x, params, Vdydz));
    dFdX->Set(5, 0, -GravPotWithParams(x, params, Vdxdz)), dFdX->Set(5, 1, -GravPotWithParams(x, params, Vdydz)), dFdX->Set(5, 2, -GravPotWithParams(x, params, Vdzdz));
    return dFdX;
}

void RightPart(double* x, double* vec, double JD, Matrix<double> *params) {
    vec[0] = x[3];
    vec[1] = x[4];
    vec[2] = x[5];
    vec[3] = 0;
    vec[4] = 0;
    vec[5] = 0;

    double rotateMatrix[3][3];

    iauC2t06a(JD + (37.0 + 32.184) / 86400.0, 0, JD, 0, 0, 0, rotateMatrix);

    for(int i = 0; i < 174 / 3; i++)
        changeCoords(rotateMatrix, x, 3 * i);
    
    double *grad = new double[3];
    grad[0] = -GravPotWithParams(x, params, Vdx);
    grad[1] = -GravPotWithParams(x, params, Vdy);
    grad[2] = -GravPotWithParams(x, params, Vdz);

    Transposition(rotateMatrix);

    changeCoords(rotateMatrix, grad, 0);

    for (int i = 0; i < 3; i++) {
        vec[i + 3] = grad[i];
    }
}

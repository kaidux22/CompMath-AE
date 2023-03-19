#include "LegFunc.h"

#define NZ_CONST 4
#define NT_CONST 4
#define NU_CONST 398600.4415 // км^3/с^2
#define R_CONST 6378.1363 // км

double GravPot(double x, double y, double z) {
	assert(y != 0);
	int Nz = NZ_CONST, Nt = NT_CONST;
	double nu = NU_CONST, R = R_CONST;
	double Jn[5] = {0.0, 0.0, -0.1082635854e-2, 0.2532435346e-5, 0.1619331205e-5};
	//первый индекс - m, второй индекс - n
	double Cmn[5][5] = { {0.0, 0.0, 0.0, 0.0, 0.0},
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

	LegFunc Pmn = LegFunc(max(Nz, Nt), Nt, sqrt(x * x + y * y) / sqrt(x * x + y * y + z * z));
	//Pmn.PrintMaxtrix();
	double res = 1;

	for (int n = 2; n <= Nz; n++) {
		res += Jn[n] * Pmn.ExtractValue(n, 0) / pow(sqrt(x * x + y * y + z * z) / R, n);
	}

	for (int n = 2; n <= Nt; n++) {
		for (int m = 1; m <= n; m++) {
			res += Pmn.ExtractValue(n, m) * (Cmn[m][n] * cos(m * atan(x / y)) + Smn[m][n]* sin(m * atan(x / y))) / pow(sqrt(x * x + y * y + z * z) / R, n);
		}
	}

	return res * (- nu / sqrt(x * x + y *y + z * z)); 
}

#define NZ_CONST 4
#define NT_CONST 4
#define NU_CONST 398600.4415 // км^3/с^2

double GravPot(double x, double y, double z) {
	int Nz = NZ_CONST, Nt = NT_CONST;
	double nu = NU_CONST;
	double Jn[5] = {0.0, 0.0, -0.1082635854e-2, 0.2532435346e-5, 0.1619331205e-5};
	double C[5][5] = { {0.0, 0.0, 0.0, 0.0, 0.0},
					   {0.0, 0.0, -0.3504890360e-9, 0.2192798802e-5, -0.5087253036e-6},
					   {0.0, 0.0, 0.1574536043e-5, 0.3090160446e-6, 0.7841223074e-7},
					   {0.0, 0.0, 0.0, 0.1005588574e-6, 0.5921574319e-7},
					   {0.0, 0.0, 0.0, 0.0, -0.3982395740e-8} };
	return 1e-5;
}

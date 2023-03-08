#include <cmath>

double* LegPol(double arg, int ord) {
	double* pols = new double[ord];
	pols[0] = 1;
	
	if (ord == 0)
		return pols;

	pols[1] = arg;

	for (int i = 1; i < ord - 1; i++) {
		pols[i + 1] = ((double)(2 * i + 1) * arg * pols[i] - (double)(i)*pols[i - 1]) / (double)(i + 1);
	}
	return pols;
}

double **LegFunc(double arg, int ord, int pwr) {
	double** funcs = new double*[pwr];

	funcs[0] = LegPol(arg, ord);


	for (int i = 1; i < pwr; i++) {
		funcs[i] = new double[ord];
		funcs[i][0] = 0;
	}

	if (ord == 0 || pwr == 0) {
		return funcs;
	}

	funcs[1][1] = sqrt(1 - arg * arg);

	for (int i = 1; i < ord - 1; i++) {

	}
	


	return funcs;
}

#include "Vnm.h"

ComplexNum Vnm(LegFunc& Pmn, int n, int m, double* vec) {
	double r = sqrt(vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2]);
	double sin = vec[2] / r;
	double cos = sqrt(vec[0] * vec[0] + vec[1] * vec[1]) / r;
	double Coef = Pmn.ExtractValue(n, m);
	Coef /= pow(r, (double)(n + m + 1));
	Coef /= pow(cos, m);
	return (ComplexNum)(Coef)*ComplexNum(vec[0], vec[1]).Pow(m);
}

ComplexNum Vdx(LegFunc& Pmn, int n, int m, double* vec) {
	if (m == 0) {
		ComplexNum V = Vnm(Pmn, n + 1, 1, vec);
		return (ComplexNum)0.5 * (-V - V.Conj());
	}
	ComplexNum res = Vnm(Pmn, n + 1, m - 1, vec) * ComplexNum((n - m + 1) * (n - m + 2));
	res = res - Vnm(Pmn, n + 1, m + 1, vec);
	return (ComplexNum)(0.5) * res;
}

ComplexNum Vdy(LegFunc& Pmn, int n, int m, double* vec) {
	if (m == 0) {
		ComplexNum V = Vnm(Pmn, n + 1, 1, vec);
		return ComplexNum(0, 0.5) * (V - V.Conj());
	}
	ComplexNum res = Vnm(Pmn, n + 1, m - 1, vec) * (ComplexNum)((n - m + 1) * (n - m + 2));
	res = res + Vnm(Pmn, n + 1, m + 1, vec);
	return ComplexNum(0, 0.5) * res;
}

ComplexNum Vdz(LegFunc& Pmn, int n, int m, double* vec) {
	return (ComplexNum)(-(n - m + 1)) * Vnm(Pmn, n + 1, m, vec);
}
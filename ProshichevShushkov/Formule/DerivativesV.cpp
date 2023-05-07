#include "DerivativesV.h"

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

ComplexNum Vdxdx(LegFunc& Pmn, int n, int m, double* vec){
    ComplexNum res;
    if(m > 1){
        res = 0.25 * Vnm(Pmn, n+2, m+2, vec);
        res -= 0.5 * (n-m+1) * (n-m+2) * Vnm(Pmn, n+2, m, vec);
        res += 0.25 * (n-m+1) * (n-m+2) * (n-m+3) * (n-m+4) * Vnm(Pmn, n+2, m-2, vec);
        return res;
    } else if(m == 1){
        res = 0.25 * Vnm(Pmn, n+2, 3, vec);
        res -= 0.5 * n * (n+1) * Vnm(Pmn, n+2, 1, vec);
        res -= 0.25 * n * (n+1) * Vnm(Pmn, n+2, 1, vec).Conj();
        return res;
    }
    res = 0.25 * Vnm(Pmn, n+2, 2, vec);
    res -= (factorial(n+2) / factorial(2*n)) * Vnm(Pmn, n+2, 0, vec);
    res += 0.25 * Vnm(Pmn, n+2, 2, vec).Conj();
    return res;
}

ComplexNum Vdydy(LegFunc& Pmn, int n, int m, double* vec){
    ComplexNum res;
    if(m > 1){
        res = -0.25 * Vnm(Pmn, n+2, m+2, vec);
        res -= 0.5 * (n-m+1) * (n-m+2) * Vnm(Pmn, n+2, m, vec);
        res -= 0.25 * (n-m+1) * (n-m+2) * (n-m+3) * (n-m+4) * Vnm(Pmn, n+2, m-2, vec);
        return res;
    } else if(m == 1){
        res = -0.25 * Vnm(Pmn, n+2, 3, vec);
        res -= 0.5 * n * (n+1) * Vnm(Pmn, n+2, 1, vec);
        res += 0.25 * n * (n+1) * Vnm(Pmn, n+2, 1, vec).Conj();
        return res;
    }
    res = -0.25 * Vnm(Pmn, n+2, 2, vec);
    res -= (factorial(n+2) / factorial(2*n)) * Vnm(Pmn, n+2, 0, vec);
    res -= 0.25 * Vnm(Pmn, n+2, 2, vec).Conj();
    return res;
}

ComplexNum Vdzdz(LegFunc& Pmn, int n, int m, double* vec){
    return (n-m+1) * (n-m+2) * Vnm(Pmn, n+2, m, vec);
}

ComplexNum Vdxdy(LegFunc& Pmn, int n, int m, double* vec){
    ComplexNum res;
    if(m > 1){
        res = Vnm(Pmn, n+2, m+2, vec) * ComplexNum(0,-0.25);
        res += ComplexNum(0,0.25) * (n-m+1) * (n-m+2) * (n-m+3) * (n-m+4) * Vnm(Pmn, n+2, m-2, vec);
        return res;
    } else if(m == 1){
        res = Vnm(Pmn, n+2, 3, vec) * ComplexNum(0,-0.25);
        res -= ComplexNum(0,0.25) * n * (n-m+1) * Vnm(Pmn, n+2, 1, vec).Conj();
        return res;
    }
    res = ComplexNum(0,-0.25) * Vnm(Pmn, n+2, 2, vec);
    res += ComplexNum(0,0.25) * Vnm(Pmn, n+2, 2, vec).Conj();
    return res;
}

ComplexNum Vdxdz(LegFunc& Pmn, int n, int m, double* vec){
    ComplexNum res;
    if(m > 0){
        res = 0.5 * (n-m+1) * Vnm(Pmn, n+2, m+1, vec);
        res -= 0.5 * (n-m+1) * (n-m+2) * (n-m+3) * Vnm(Pmn, n+2, m-1, vec);
        return res;
    }
    res = 0.5 * (n+1) * Vnm(Pmn, n+2, 1, vec);
    res += 0.5 * (n+1) * Vnm(Pmn, n+2, 1, vec).Conj();
    return res;
}

ComplexNum Vdydz(LegFunc& Pmn, int n, int m, double* vec){
    ComplexNum res;
    if(m > 0){
        res = ComplexNum(0,-0.5) * (n-m+1) * Vnm(Pmn, n+2, m+1, vec);
        res -= ComplexNum(0,0.5) * (n-m+1) * (n-m+2) * (n-m+3) * Vnm(Pmn, n+2, m-1, vec);
        return res;
    }
    res = ComplexNum(0,-0.5) * (n+1) * Vnm(Pmn, n+2, 1, vec);
    res += ComplexNum(0,0.5) * (n+1) * Vnm(Pmn, n+2, 1, vec).Conj();
    return res;
}


#ifndef COMPMATH_VNM
#define COMPMATH_VNM

#include "ComplexNums.h"
#include "LegFunc.h"

//функция для дальнейшего подсчёта градиента реккурентно
ComplexNum Vnm(LegFunc& Pmn, int n, int m, double x, double y, double z) {
	double Coef = Pmn.ExtractValue(n, m) / (pow(sqrt(x * x + y * y + z * z), (double)(n + m + 1)));
	return (ComplexNum)(Coef)*ComplexNum(x, y).Pow(m);
}

ComplexNum Vdx(LegFunc& Pmn, int n, int m, double x, double y, double z) {
	return (ComplexNum)((-Vnm(Pmn, n + 1, m + 1, x, y, z)).Real());
}

ComplexNum Vdy(LegFunc& Pmn, int n, int m, double x, double y, double z) {
	return (ComplexNum)((-Vnm(Pmn, n + 1, m + 1, x, y, z)).Img());
}

ComplexNum Vdz(LegFunc& Pmn, int n, int m, double x, double y, double z) {
	return (ComplexNum)(- (n - m + 1)) * Vnm(Pmn, n + 1, m, x, y, z);
}

#endif //COMPMATH_VNM
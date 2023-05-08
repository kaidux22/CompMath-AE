#ifndef COMPMATH_INTEGRATEFUNC
#define COMPMATH_INTEGRATEFUNC

#include "../Formule/GravPot.h"

/* Численное интегрирование методом Дорманда-Принца */
void DormandPrince(double JD, double h, const int N, double* vec, double a[7][7], double b[7], double** k, double c[7], void (*f)(double*, double*, double));

/* Решение систему ОДУ */
double **intergrate(double JD_start, double h, const int N, double *vec);

#endif //COMPMATH_INTEGRATEFUNC
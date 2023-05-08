#ifndef COMPMATH_RIGHTPART
#define COMPMATH_RIGHTPART

#include "../Math/Matrix.h"
#include "../sofa/sofa.h"
#include "../Formule/Converter.h"
#include "../Formule/GravPot.h"

/* Гравитационный потенцил с кастомизацией параметров */
double GravPotWithParams(double* vec, Matrix<double> *params, ComplexNum(*func)(LegFunc&, int, int, double*));

/* Построение матрицы dF/dX для численного интегрирования */
Matrix<double> *MatrixdFdX(double *x, Matrix<double> *params);

/* Правая часть ОДУ второго порядка */
void RightPart(double* x, double* vec, double JD, Matrix<double> *params);

#endif //COMPMATH_RIGHTPART
#ifndef COMPMATH_RIGHTPART
#define COMPMATH_RIGHTPART

#include "../sofa/sofa.h"
#include "../Formule/Converter.h"
#include "../Formule/GravPot.h"

/* Построение матрицы df/dx для численного интегрирования */
Matrix<double> *MatrixdFdX(double *x, Matrix<double> *params, double JD);

/* Построение матрицы df/dp для численного интегрирования */
Matrix<double> *MatrixdFdParam(double *x, Matrix<double> *params, double JD);

/* Правая часть ОДУ второго порядка в обратной задаче */
void RightPart(double* x, double* vec, double JD, Matrix<double> *params);

#endif //COMPMATH_RIGHTPART
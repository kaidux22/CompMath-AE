#ifndef COMPMATH_GRAVPOT
#define COMPMATH_GRAVPOT

#include <cmath>
#include <cassert>
#include <iostream>

#include "../Math/LegFunc.h"
#include "../Math/ComplexNums.h"
#include "Vnm.h"
#include "Converter.h"
#include "../sofa/sofa.h"

#define N_CONST 4
#define NU_CONST 398600.4415 // км^3/с^2
#define R_CONST 6378.1363 // км
#define MAX_ORD 2 //наибольшая степень производной
#define GENERAL_TIME 86400.0 //сутки в секундах

using namespace std;

/* Вычисление гравитационного потенциала */
double GravPot(double* vec, ComplexNum(*func)(LegFunc&, int, int, double*));

/* Производная по Cnm */
double DerivativedVdC(double* vec, int n, int m, ComplexNum(*func)(LegFunc&, int, int, double*));

/* Производная по Snm */
double DerivativedVdS(double* vec, int n, int m, ComplexNum(*func)(LegFunc&, int, int, double*));

/* Производная по M */
double DerivativedVdM(double* vec, ComplexNum(*func)(LegFunc&, int, int, double*));

/* Подсчёт градиента гравитационного потенциала */
void GradV(double* x, double* vec, double JD);

/* Численное интегрирование методов Дорманда-Принца */
void DormandPrince(double JD, double h, const int N, double* vec, double a[7][7], double b[7], double** k, double c[7], void (*f)(double*, double*, double));

/* Решение систему ОДУ */
double **intergrate(double JD_start, double h, const int N, double* vec);

#endif //COMPMATH_GRAVPOT
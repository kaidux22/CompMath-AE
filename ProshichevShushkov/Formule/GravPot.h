#ifndef COMPMATH_GRAVPOT
#define COMPMATH_GRAVPOT

#include <cmath>
#include <cassert>
#include <iostream>

#include "../Math/LegFunc.h"
#include "../Math/ComplexNums.h"
#include "DerivativesV.h"
#include "Converter.h"
#include "../sofa/sofa.h"

#define N_CONST 2
#define NU_CONST 398600.4415 // км^3/с^2
#define R_CONST 6378.1363 // км
#define M_CONST 5,972e24 // кг
#define MAX_ORD 2 //наибольшая степень производной

using namespace std;

/* Вычисление гравитационного потенциала */
double GravPot(double* vec, ComplexNum(*func)(LegFunc&, int, int, double*));

/* Подсчёт градиента гравитационного потенциала */
double* GradV(double* vec, double UTC);

/* Численное интегрирование методов Дорманда-Принца */
void DormandPrince(double UTC, double h, const int N, double* vec, double a[7][7], double b[7], double** k, int integrate_numder, double* (*f)(double* vec, double));

/* Решение систему ОДУ */
double **intergrate(double UTC_start, double h, const int N, double* vec);

double DerivativedVdC(double* vec, ComplexNum(*func)(LegFunc&, int, int, double*));

double DerivativedVdS(double* vec, int n, int m, ComplexNum(*func)(LegFunc&, int, int, double*));

#endif //COMPMATH_GRAVPOT
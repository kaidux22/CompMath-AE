#ifndef COMPMATH_GRAVPOT
#define COMPMATH_GRAVPOT

#include <cmath>
#include <cassert>
#include <iostream>

#include "LegFunc.h"
#include "ComplexNums.h"
#include "Vnm.h"
#include "Converter.h"

#define N_CONST 2
#define NU_CONST 398600.4415 // км^3/с^2
#define R_CONST 6378.1363 // км
#define MAX_ORD 2 //наибольшая степень производной

using namespace std;

/* Вычисление гравитационного потенциала */
double GravPot(double* vec, ComplexNum(*func)(LegFunc&, int, int, double*));

/* Подсчёт градиента гравитационного потенциала */
double* GradV(double* vec, double UTC);

/* Численное интегрирование методов Дорманда-Принца */
void DormandPrince(double t, double h, const int N, double* vec, double a[7][7], double b[7], int integrate_numder, double* (*f)(double* vec, double));

/* Решение систему ОДУ */
void intergrate(double UTC_start, double h, const int N, double* vec);

#endif //COMPMATH_GRAVPOT
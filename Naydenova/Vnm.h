#pragma once
#include "ComplexNums.h"
#include "LegFunc.h"

//функция для дальнейшего подсчёта градиента реккурентно
ComplexNum Vnm(LegFunc& Pmn, int n, int m, double *vec) ;

ComplexNum Vdx(LegFunc& Pmn, int n, int m, double* vec) ;

ComplexNum Vdy(LegFunc& Pmn, int n, int m, double* vec) ;

ComplexNum Vdz(LegFunc& Pmn, int n, int m, double* vec) ;
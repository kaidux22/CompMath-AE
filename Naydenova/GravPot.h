#pragma once
#include <cmath>
#include <cassert>
#include <iostream>
#include "Vnm.h"
#include "sofa/sofa.h"
#include "Converter.h"

#define N_CONST 2
#define NU_CONST 398600.4415 // км^3/с^2
#define R_CONST 6378.1363 // км
#define MAX_ORD 2

using namespace std;

double GravPot(double *vec, ComplexNum(*d)(LegFunc&, int, int, double*)) ;

void GradV(double* vec, double time) ;

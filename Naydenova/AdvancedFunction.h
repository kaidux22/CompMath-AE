#pragma once

#include <cmath>
#include <cassert>
#include <iostream>
#include "DoubleDerivatives.h"
#include "Derivatives.h"
#include "sofa/sofa.h"
#include "Converter.h"

using namespace std;

double** create_matrix_df_dx(double* x, double mu, double J, double JD);
void function(double* x, double* vec, double time, double J, double mu) ;
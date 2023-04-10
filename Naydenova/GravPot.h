#pragma once
#include <cmath>
#include <cassert>
#include <iostream>
#include "Derivatives.h"
#include "sofa/sofa.h"
#include "Converter.h"

#define R_CONST 6378.1363 // км

using namespace std;

void GradV(double* vec, double time) ;

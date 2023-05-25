#pragma once

#include <iostream>
#include <cmath>
#include <cassert>
#include "Observatories.h"
#include "sofa/sofa.h"
#include "DormandPrince.h"
#include "GravPot.h"
#include "Converter.h"
#include "sofa/sofa.h"
#include <fstream>
#include "SLE.h"
#include "Derivatives.h"
#include <random>
#include <vector>

using namespace std;

double** create_observatories(double JD_start);
double** create_model_observatories(double JD_start);
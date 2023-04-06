#pragma once

#include <iostream>
#include <cmath>
#include <cassert>

#include "Formule/GravPot.h"
#include "Math/ComplexNums.h"
#include "Formule/Converter.h"
#include "sofa/sofa.h"

/* Юлианская дата */
#define JD 2451545.0

/* Шаг в 1 минуту */
#define STEP 30.0 //с

/* G * m */
#define GM 398600.4415 // км^3 / с^2

/* радиус Земли + 500 км */
#define START_POINT 6878.0 //км




#pragma once
#include <cmath>
#include <iostream>
using namespace std;
#define PI 3.1415926535

#define R_CONST 6378.1363 // км

class Observatories {
public:
    Observatories(double longitude, double cos, double sin);
    double* get_coords();

private:
    void SphericalCoordToCartesian();
    double longitude;
    double angle;
    double cos_l;
    double sin_l;
    double* vec;
};
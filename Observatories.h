#pragma once
#include <cmath>
#include <iostream>
using namespace std;

#define R_CONST 6378.1363 // км

class Observatories {
public:
	Observatories(double longitude, double cos, double sin): longitude(longitude), cos_l(cos), sin_l(sin) {
		vec = new double[3];
	};
	void СylindrСoordToCartesian();
	double* get_coords();

private:
	double longitude;
	double cos_l;
	double sin_l;
	double* vec;
};

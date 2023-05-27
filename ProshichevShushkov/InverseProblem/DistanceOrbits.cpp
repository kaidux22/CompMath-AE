#include "DistanceOrbits.h"
#include <iostream>

double* OrbitDistance(double **orbits, int stepCnt){
	double *distance = new double[stepCnt * 2];
	for(int i = 0; i < stepCnt; i++){
		distance[2 * i] = orbits[i][0];
		distance[2 * i + 1] =  sqrt(pow(orbits[i][1] - orbits[i][4], 2) +
								   pow(orbits[i][2] - orbits[i][5], 2) +
								   pow(orbits[i][3] - orbits[i][6], 2));
	}
	return distance;
}
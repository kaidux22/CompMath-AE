#include "DistanceOrbits.h"
#include <iostream>

double** OrbitDistance(double** orbit1, double** orbit2, int stepCnt){
	double** distance = new double*[stepCnt];
	for(int i = 0; i < stepCnt; i++){
		distance[i] = new double[2];
		distance[i][0] = orbit1[i][0];
		distance[i][1] = sqrt(pow(orbit1[i][1] - orbit2[i][1], 2) +
				   pow(orbit1[i][2] - orbit2[i][2], 2) +
				   pow(orbit1[i][3] - orbit2[i][3], 2));
	}
	return distance;
}
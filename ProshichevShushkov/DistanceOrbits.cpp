#include "DistanceOrbits.h"
#include <iostream>

double** OrbitDistance(double** orbit1, double** orbit2, int stepCnt){
	double* distance = new double[stepCnt];
	for(int i = 1; i < stepCnt; i++)
		distance[i] = sqrt(pow(orbit1[i][1] - orbit2[i][1], 2) +
				   pow(orbit1[i][2] - orbit2[i][2], 2) +
				   pow(orbit1[i][3] - orbit2[i][3], 2));

	double** distancePerTime = new double*[stepCnt + 1];
	for(int i = 1; i < stepCnt - 1; i++){
		distancePerTime[i - 1] = new double[2];
		distancePerTime[i - 1][0] = orbit1[i][0];
		distancePerTime[i - 1][1] = (distance[i] - distance[i - 1]) / (orbit1[i][0] - orbit1[i - 1][0]);
	}	
	delete[] distance;
	return distancePerTime;
}
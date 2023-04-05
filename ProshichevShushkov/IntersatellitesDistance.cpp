#include <iostream>
#include "GravPot.h"
#include "Converter.h"

#define GENERAL_TIME 86400

void ChangeCoordsTypeToEarth(double** vec, int cnt, double JD){
    double rotateMatrix[3][3];
    iauC2t06a(JD + (37.0 + 32.184) / 86400.0, 0, JD, 0, 0, 0, rotateMatrix);
    for(int i = 0; i < cnt; ++i){
        double* newCoords = new double[4];
        newCoords[0] = vec[i][0];
        for (int j = 0; j < 3; j++) {
            newCoords[j+1] = rotateMatrix[j][0] * vec[i][1]
                    + rotateMatrix[j][1] * vec[i][2]
                    + rotateMatrix[j][2] * vec[i][3];
        }
        for(int j = 0; j < 3; ++j)
            vec[i][j+1] = newCoords[i+1];
        delete newCoords;
    }
}


float** IntersatelliteDistanceChange(double** FirstSatelliteOrbit, double** SecondSatelliteOrbit, double JD, double h){
    int cnt = GENERAL_TIME / h;

    ChangeCoordsTypeToEarth(FirstSatelliteOrbit, cnt, JD);
    ChangeCoordsTypeToEarth(SecondSatelliteOrbit, cnt, JD);

    float* IntersatelliteDistance = new float[cnt];
    for(int i = 1; i < cnt; ++i){
        IntersatelliteDistance[i] = sqrt(pow(FirstSatelliteOrbit[i][1] - SecondSatelliteOrbit[i-1][1], 2.0)
                +pow(FirstSatelliteOrbit[i][2] - SecondSatelliteOrbit[i-1][2], 2.0)
                +pow(FirstSatelliteOrbit[i][3] - SecondSatelliteOrbit[i-1][3], 2.0));
    }

    float** DistanceChangePerTime = new float*[cnt-1];
    for(int i = 1; i < cnt; ++i){
        DistanceChangePerTime[i-1] = new float[2];
        DistanceChangePerTime[i-1][0] = FirstSatelliteOrbit[i][0];
        DistanceChangePerTime[i-1][0] = (IntersatelliteDistance[i] - IntersatelliteDistance[i-1])
                /(FirstSatelliteOrbit[i][0]-FirstSatelliteOrbit[i-1][0]);
    }

    delete[] IntersatelliteDistance;
    return DistanceChangePerTime;
}
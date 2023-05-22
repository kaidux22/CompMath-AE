#include "Converter.h"

void Transposition(double(*matrix)[3]) {
    for (int i = 0; i < 3; i++) {
        for (int j = i; j < 3; j++) {
            double tmp = matrix[i][j];
            matrix[i][j] = matrix[j][i];
            matrix[j][i] = tmp;
        }
    }
}

void changeCoords(double(*rotateMatrix)[3], double* vec, int idx) {
    double* newCoords = new double[3];
    for (int i = 0; i < 3; i++) {
        newCoords[i] = rotateMatrix[i][0] * vec[idx] + rotateMatrix[i][1] * vec[idx + 1] + rotateMatrix[i][2] * vec[idx + 2];
    }

    for(int i = 0; i < 3; i++)
        vec[idx + i] = newCoords[i];
    delete newCoords;
}

void changeCoordsRight(double(*rotateMatrix)[3], double* vec, int idx) {
    double* newCoords = new double[3];
    for (int i = 0; i < 3; i++) {
        newCoords[i] = rotateMatrix[0][i] * vec[idx] + rotateMatrix[1][i] * vec[idx + 1] + rotateMatrix[2][i] * vec[idx + 2];
    }

    for(int i = 0; i < 3; i++)
        vec[idx + i] = newCoords[i];
    delete newCoords;
}

double** multiplication_AtA(vector<vector<double>> &A){
    double** AtA = new double * [8];
    for (int i=0 ; i < 8; i++){
        AtA[i] = new double [8];
    }

    int n = A.size();

    for(int i=0; i < 8; i++){
        for (int j=0; j < 8; j++){
            double res = 0;
            for (int t = 0; t < n; t++){
                res += A[t][i] * A[t][j];
            }
            AtA[i][j] = res;
        }
    }
    return AtA;
}

double* multiplication_Atr(vector<vector<double>> &A, vector<double>& r){
    double* Atr = new double[8];
    int n = A.size();

    for(int i=0; i < 8; i++){
        double res = 0;
        for (int t = 0; t < n; t++){
            res += A[t][i] * r[t];
        }
        Atr[i] = res;
    }
    return Atr;
}
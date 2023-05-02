#include "SLE.h"

// Ax = b
// LL^t * x = b
// L^t * x = y
// Ly = b

double* Cholesky_decomposition(double** A, int size, double* b){
    double **L = new double* [size];
    for (int i=0; i < size; i++){
        L[i] = new double [size];
    }

    for (int i=0; i < size; i++){
        for (int j = 0; j < (i + 1); j++){
            double sum = 0;
            for (int k = 0; k < j; k++) {
                sum += L[i][k] * L[j][k];
            }

            if (i == j) {
                L[i][j] = sqrt(A[i][i] - sum);
            } else {
                L[i][j] = (1.0 / L[j][j] * (A[i][j] - sum));
            }

        }
    }
    double* x = new double[size];

    //  L*y=b
    double* y = new double[size];
    for (int i = 0; i < size; i++){
        double sum = 0;
        for (int j = 0; j < i; j++){
            sum += L[i][j] * y[j];
        }
        y[i] = (1.0 / L[i][i]) * (b[i] - sum);
    }

    //  L^t*x=y
    for (int i = size-1; i >= 0; i--){
        double sum = 0;
        for (int j = i+1; j < size; j++){
            sum += L[j][i] * x[j];
        }
        x[i] = (1.0 / L[i][i]) * (y[i] - sum);
    }

    delete[] y;

    for (int i=0; i < size; i++){
        delete[] L[i];
    }
    delete[] L;

    return x;
}


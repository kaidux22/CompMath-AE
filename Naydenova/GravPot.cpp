#include "GravPot.h"



void GradV(double* x, double* vec, double JD, double J, double mu) {

    vec[0] = x[3];
    vec[1] = x[4];
    vec[2] = x[5];
    vec[3] = 0;
    vec[4] = 0;
    vec[5] = 0;

    double rotateMatrix[3][3];
    /*
    for(int i=0; i < 3; i++){
        for (int j =0; j < 3; j++){
            rotateMatrix[i][j] = 0;
        }
    }

    rotateMatrix[0][0] = 1;
    rotateMatrix[1][1] = 1;
    rotateMatrix[2][2] = 1; */
    iauC2t06a(JD + (37.0 + 32.184) / 86400.0, 0, JD, 0, 0, 0, rotateMatrix);

    changeCoords(rotateMatrix, x, 0);   // НСК -> ЗСК

    double *grad = new double[3];
    grad[0] = -dx(x, MU_CONST, J2);
    grad[1] = -dy(x, MU_CONST, J2);
    grad[2] = -dz(x, MU_CONST, J2);

    Transposition(rotateMatrix);

    changeCoords(rotateMatrix, grad, 0); // ЗСК -> НСК

    for (int i = 0; i < 3; i++) {
        vec[i + 3] = grad[i];
    }
}


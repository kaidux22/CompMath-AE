#include "GravPot.h"



void GradV(double* vec, double JD) {

    double rotateMatrix[3][3];

    iauC2t06a(JD + (37.0 + 32.184) / 86400.0, 0, JD, 0, 0, 0, rotateMatrix);

    changeCoords(rotateMatrix, vec, 0);

    double* grad = new double[3];
    grad[0] = dx(vec);
    grad[1] = dy(vec);
    grad[2] = dz(vec);

    Transposition(rotateMatrix);

    changeCoords(rotateMatrix, grad, 0);

    for (int i=0; i < 3; i++) {
        vec[i] = vec[i+3];
        vec[i+3] = grad[i];
    }

}

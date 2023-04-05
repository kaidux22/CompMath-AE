#include "GravPot.h"



double GravPot(double *vec, ComplexNum(*d)(LegFunc&, int, int, double*)) {
    assert(vec[1] != 0);

    int N = N_CONST;
    double nu = NU_CONST, R = R_CONST;

    double J0 = 1.0; double J2 = -0.1082635854e-2;

    LegFunc Pmn = LegFunc(N + MAX_ORD, N + MAX_ORD, sqrt(vec[0] * vec[0] + vec[1] * vec[1]) / sqrt(vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2]));
    ComplexNum res = ComplexNum(0, 0);

    res = (ComplexNum)J0 * d(Pmn, 0, 0, vec) + (ComplexNum)pow(R, 2) * (ComplexNum)J2 * d(Pmn, 2, 0, vec);

    return res.Real() * (-nu);

}

void GradV(double* vec, double JD) {

    double rotateMatrix[3][3];

    iauC2t06a(JD + (37.0 + 32.184) / 86400.0, 0, JD, 0, 0, 0, rotateMatrix);

    changeCoords(rotateMatrix, vec);

    double* grad = new double[6];
    grad[0] = GravPot(vec, Vdx);
    grad[1] = GravPot(vec, Vdy);
    grad[2] = GravPot(vec, Vdz);

    Transposition(rotateMatrix);

    changeCoords(rotateMatrix, grad);

    for (int i=0; i < 3; i++) {
        vec[i] = vec[i+3];
        vec[i+3] = grad[i];
    }

}

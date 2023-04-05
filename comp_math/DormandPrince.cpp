#include "DormandPrince.h"
/*
y'(t) = v(t)
v'(t) = 2
 */

void DormandPrince(double JD, double h, const int N, double* vec, double a[7][7], double b[7], double** k, void (*f)(double*, double)) {

    for (int i = 0 ; i < 7; i++) {
        for (int j = 0; j < N; j++) {
            k[i][j] = vec[j];
            for (int t = 0; t < 7; t++) {
                k[i][j] += a[i][t] * h * k[t][j];
            }
        }
        f(k[i], JD);
    }

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < 7; j++) {
            vec[i] += b[j] * k[j][i] * h;
        }
    }

}


double** intergrate(double JD, double h, const int N, double* vec) {

    double a[7][7] = { {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
                       {1.0 / 5.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
                       {3.0 / 40.0, 9.0 / 40.0, 0.0, 0.0, 0.0, 0.0, 0.0},
                       {44.0 / 45.0, -56.0 / 15.0, 32.0 / 9.0, 0.0, 0.0, 0.0, 0.0 },
                       {19372.0 / 6561.0, -25360.0 / 2187.0, 64448.0 / 6561.0, -212.0 / 729.0, 0.0, 0.0, 0.0 },
                       {9017.0 / 3168.0, -355.0 / 33.0, 46732.0 / 5247.0, 49.0 / 176.0, -5103.0 / 18656.0, 0.0, 0.0},
                       {35.0 / 384.0, 0.0, 500.0 / 1113.0, 125.0 / 192.0, -2187.0 / 6784.0, 11.0 / 84.0, 0.0},
    };

    double b[7] = { 5179.0 / 57600.0, 0.0, 7571.0 / 16695.0, 393.0 / 640.0, -92097.0 / 339200.0, 187.0 / 2100.0, 1.0 / 40.0 };

    double** k = new double* [7];
    for (int i = 0; i < 7; i++) {
        k[i] = new double[N];
        for (int j = 0; j < N; j++) {
            k[i][j] = 0;
        }
    }

    int cnt = GENERAL_TIME / h;
    double** orbit = new double* [cnt];

    // 86400 секунд в сутках
    for (int i = 0; i < cnt; i++) {
        DormandPrince(JD, h, N, vec, a, b, k, GradV);
        orbit[i] = new double[7];
        orbit[i][0] = JD, orbit[i][1] = vec[0], orbit[i][2] = vec[1], orbit[i][3] = vec[2], orbit[i][4] = vec[3], orbit[i][5] = vec[4], orbit[i][6] = vec[5];
        JD += h / 86400.0;
    }

    for (int i = 0; i < N; i++) {
        delete[] k[i];
    }
    delete[] k;

    return orbit;

}



/*
int main() {

    int N = 3;
    auto* y = new double[N];
    y[0] = 1.0; y[1] = 0.5; y[2] = 0.1;
    double t = 0.0, h = 0.01;
    intergrate(t, h, N, y);


    return 0;
}
*/
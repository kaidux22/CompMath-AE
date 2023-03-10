#include <iostream>
#include <cmath>

using namespace std;


float func(float x){
    return -5*x;
}


long double DormanPrince(long double x0, long double y0, long double h, int x){
    long double a[7][6] = {{0, 0, 0, 0, 0, 0},
             {1 / 5, 0, 0, 0, 0, 0},
             {3 / 40, 9 / 40, 0, 0, 0, 0},
             {44 / 45, -56 / 15, 32 / 9, 0, 0, 0},
             {19372 / 6561, -25360 / 2187, 64448 / 6561, -212 / 729, 0, 0},
             {9017 / 3168, -355 / 33, 46732 / 5247, 49 / 176, -5103 / 18656, 0},
             {35 / 384, 0, 500 / 1113, 125 / 192, -2187 / 6784, 11 / 84}};

    long double b1[7] = {5179 / 57600, 0, 7571 / 16695, 393 / 640, -92097 / 339200, 187 / 2100, 1 / 40};
    long double b[7] = {35 / 384, 0, 500 / 1113, 125 / 192, -2187 / 6784, 11 / 84, 0};
    long double k[7];

    long double yk = y0, zk = y0;
    long double eps = 0.000001, err, s, hopt;
    int i = 0;

    while (x0 < x){
        k[0] = h * func(yk);
        k[1] = h * func( yk + a[1][0]* k[0]);
        k[2] = h * func( yk + a[2][0] * k[0] + a[2][1] * k[1]);
        k[3] = h * func( yk + a[3][0] * k[0] + a[3][1] * k[1] + a[3][2] * k[2]);
        k[4] = h * func( yk + a[4][0] * k[0] + a[4][1] * k[1] + a[4][2] * k[2] + a[4][3] * k[3]);
        k[5] = h * func( yk + a[5][0] * k[0] + a[5][1] * k[1] + a[5][2] * k[2] + a[5][3] * k[3] + a[5][4] * k[4]);
        k[6] = h * func( yk + a[6][0] * k[0] + a[6][1] * k[1] + a[6][2] * k[2] + a[6][3] * k[3] + a[6][4] * k[4] + a[6][5] * k[5]);

        for (int j =0; j < 7; j++){
            yk += k[j] *b[j];
            zk += k[j] *b1[j];
        }

        err = abs(yk - zk);
        s = pow(eps*h/(2*err), 1/5);
        hopt = s * h;
        if (h <= hopt){
            x0 = x0 + h;
            y0 = yk;
            h = hopt;
        }
        else {
            h = hopt;
        }

    }
    return yk;
}


int main()
{
    float x0 = 0, y = 0, h = 0.1, x = 5;
    return 0;
}



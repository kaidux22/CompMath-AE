#include <iostream>
#include <cmath>
#include <vector>

using namespace std;


float func(float x){
    return -x;
}


std::vector<long double> DormanPrince(long double x0, long double y0, long double h, long double x){

    vector<vector <long double>> a = {{0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
                             {1.0 / 5.0, 0.0, 0.0, 0.0, 0.0, 0.0},
                             {3.0 / 40.0, 9.0 / 40.0, 0.0, 0.0, 0.0, 0.0},
                             {44.0 / 45.0, -56.0 / 15.0, 32.0 / 9.0, 0.0, 0.0, 0.0},
                             {19372.0 / 6561.0, -25360.0 / 2187.0, 64448.0 / 6561.0, -212.0 / 729.0, 0.0, 0.0},
                             {9017.0 / 3168.0, -355.0 / 33.0, 46732.0 / 5247.0, 49.0 / 176.0, -5103.0 / 18656.0, 0.0},
                             {35.0 / 384.0, 0.0, 500.0 / 1113.0, 125.0 / 192.0, -2187.0 / 6784.0, 11.0 / 84.0}};

    vector <long double> b1 = {5179.0 / 57600.0, 0.0, 7571.0 / 16695.0, 393.0 / 640.0, -92097.0 / 339200.0, 187.0 / 2100.0, 1.0 / 40.0};
    vector <long double> b = {35.0 / 384.0, 0.0, 500.0 / 1113.0, 125.0 / 192.0, -2187.0 / 6784.0, 11.0 / 84.0, 0.0};
    vector <long double> k = {0,0,0,0,0,0,0,0};

    long double yk = y0, zk = y0;
    long double tol = 0.0001, err;
    vector<long double> y; y.push_back(y0);
    long double hmax = 0.2; long double facmax = 1.0;
    long double fac = 0.8; long double facmin = 0.6;
    long double gamma;              //??? надо выбрать константы

    while (x0 < x){
        k[0] = h * func(y0);
        k[1] = h * func( y0 + a[1][0]* k[0]);
        k[2] = h * func( y0 + a[2][0] * k[0] + a[2][1] * k[1]);
        k[3] = h * func( y0 + a[3][0] * k[0] + a[3][1] * k[1] + a[3][2] * k[2]);
        k[4] = h * func( y0 + a[4][0] * k[0] + a[4][1] * k[1] + a[4][2] * k[2] + a[4][3] * k[3]);
        k[5] = h * func( y0 + a[5][0] * k[0] + a[5][1] * k[1] + a[5][2] * k[2] + a[5][3] * k[3] + a[5][4] * k[4]);
        k[6] = h * func( y0 + a[6][0] * k[0] + a[6][1] * k[1] + a[6][2] * k[2] + a[6][3] * k[3] + a[6][4] * k[4] + a[6][5] * k[5]);

        for (int j =0; j < 7; j++){ //подсчет ук и zк на очередном шаге
            yk += k[j] *b[j];
            //cout << yk;
            zk += k[j] *b1[j];
        }

        err = (yk-zk)/(sqrt(zk*zk)); //пока так(
        gamma = min(facmax, max(facmin, pow(fac*(tol/err), 1/5)));

        if (err <= tol){
            x0 = x0 + h;
            y0 = yk;
            h = min(gamma*h, hmax);
            y.push_back(yk);
        }
        else{
            h = gamma*h;
        }

    }
    return y;
}


int main(){
    long double x0 = 0, y0 = 1, h = 0.2, x = 1;
    std::vector <long double> y = DormanPrince(x0, y0, h, x);
    for (int i=0; i < y.size(); i++){
        cout << y[i] << ' ';
    }
    return 0;
}



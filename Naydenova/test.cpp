#include <vector>
#include <iostream>

using namespace std;



int main(){

    double *states = new double[48];

    double x[54];
    for (int i=0; i< 9; i++){
        for (int j=0; j < 6; j++){
            x[i*6 + j ] = i;
        }
    }

    for (int i =0; i < 54; i++){
        cout << x[i] << " ";
    }
    cout << endl;

    double df_dx[6][6] = {{0,0,0,1,0,0},
                          {0,0,0,0,1,0},
                          {0,0,0,0,0,1},
                          {1, 2, 3,0,0,0},
                          {4, 5, 6,0,0,0},
                          {7, 8, 9,0,0,0}};

    for (int i=0; i < 6; i++){
        for (int j = 0; j < 8; j++){
            double res = 0;
            for (int k = 0; k < 6; k++){
                res += df_dx[i][k] * x[6 + 6*j + k];
            }
            states[6*j + i] = res;  //НСК
        }
    }

    for (int i =0; i < 48; i++){
        cout << states[i] << " ";
    }
}
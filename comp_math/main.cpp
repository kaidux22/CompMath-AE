#include "main.h"

#define JD 2451545.0
#define STEP 60.0

int main(){

    Observatories first = Observatories(2.33675, 0.659470, 0.749223);
    Observatories second = Observatories(359.1083,  0.70862,  0.70323);
    Observatories third = Observatories(28.9667, 0.75566, 0.65278);
    double* station1 = first.get_coords();
    double* station2 = second.get_coords();
    double* station3 = third.get_coords();


    double JD_start = JD;
    double *vec = new double[6];
    double rotateMatrix[3][3];
    vec[0] = -4661.36, vec[1] = -3953.12, vec[2] = 3154.59, vec[3] = 0, vec[4] = 0, vec[5] = 0;
    iauC2t06a(JD_start + (37.0 + 32.184) / 86400.0, 0, JD_start, 0, 0, 0, rotateMatrix);
    Transposition(rotateMatrix);
    changeCoords(rotateMatrix, vec);

    double** res = intergrate(JD, STEP, 6, vec);

    for(int i = 0; i < 10; i++){
        cout << "time: " << res[i][0] << " x: " << res[i][1] << " y: " << res[i][2] << " z: " << res[i][3] ;
        cout  << " " << res[i][4] << " " << res[i][5] << " " << res[i][6] << endl;
        cout << endl;
    }
    delete[] vec;

}
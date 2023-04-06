#include "main.h"

#define JD 2451545.0
#define STEP 60.0
#define GM 398600.4415 // км^3 / с^2
#define START_POINT 6878.0 //км


double** create_observatories(){

    Observatories initial_coord[8] = {Observatories(2.33675, 0.659470, 0.749223),
                                      Observatories(359.1083,  0.70862,  0.70323),
                                      Observatories(28.9667, 0.75566, 0.65278),
                                      Observatories(260.8053,  0.94388,  0.33026),
                                      Observatories(9.10031, 0.698332, 0.713430),
                                      Observatories(149.0642,  0.85563,  0.51621),
                                      Observatories(117.5750,  0.76278,  0.64470),
                                      Observatories(342.11094, 0.877701, 0.478380),
                                      };

    double** station = new double *[8];
    for (int i=0; i < 8; i++){
        station[i] = initial_coord[i].get_coords();
    }

    double rotateMatrix[3][3];
    double JD_start = JD;
    iauC2t06a(JD_start + (37.0 + 32.184) / 86400.0, 0, JD_start, 0, 0, 0, rotateMatrix);
    Transposition(rotateMatrix);
    for (int i=0; i < 8; i++){
        changeCoords(rotateMatrix, station[i], 0);
    }

    return station;

}

int main(){

    int cnt = GENERAL_TIME / STEP;
    double JD_start = JD;
    double *vec = new double[6];
    double rotateMatrix[3][3];
    vec[0] = START_POINT, vec[1] = 0, vec[2] = 0;
    vec[3] = 0, vec[4] = sqrt(GM / START_POINT), vec[5] = 0;

    iauC2t06a(JD_start + (37.0 + 32.184) / 86400.0, 0, JD_start, 0, 0, 0, rotateMatrix);
    Transposition(rotateMatrix);
    changeCoords(rotateMatrix, vec, 0);
    changeCoords(rotateMatrix, vec, 3);

    double** res = intergrate(JD, STEP, 6, vec);

    double** stations = create_observatories();


    for (int i=0; i < cnt; i++) {
        for (int j = 0; j < 8; j++) {
            cout << sqrt(pow((stations[j][0] - res[i][1]), 2) + pow((stations[j][1] - res[i][2]), 2) +
                         pow((stations[j][2] - res[i][3]), 2)) << "    ";
        }
        cout << endl;
    }

    for(int i=0; i < cnt; i++){
        delete[] res[i];
    }
    delete[] res;

    for (int i=0; i < 8; i++){
        delete[] stations[i];
    }
    delete[] stations;

    delete[] vec;
}
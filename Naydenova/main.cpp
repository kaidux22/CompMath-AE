#include "main.h"

#define JD 2451545.0
#define STEP 30.0
#define GM 398600.4415 // км^3 / с^2
#define START_POINT 6878.0 //км


double** create_observatories(double JD_start){

    Observatories initial_coord[8] = {
                                      Observatories(281.5075,  1.00045,  -0.00405),
                                      Observatories(260.8053,  0.94388,  +0.33026),
                                      Observatories(42.5008,  0.72958,  0.68232 ),
                                      Observatories(281.65,    0.999,    +0.000),
                                      Observatories(321.3126,  0.98840,  -0.15179),
                                      Observatories(299.99039, 0.998647, -0.051941),
                                      Observatories(101.43942, 0.998617, +0.052565),
                                      Observatories(282.70896, 1.000183, +0.021030),
                                      };

    double** station = new double *[8];
    for (int i=0; i < 8; i++){
        station[i] = initial_coord[i].get_coords();
    }


    double rotateMatrix[3][3];
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

    double** res = integrate(JD, STEP, 6, vec);

    /*
    for (int i=0; i < cnt; i++) {
        cout << res[i][0] << "    " << res[i][1] << "    " << res[i][2] << "    " << res[i][3] << "    " << res[i][4] << "    " << res[i][5] << "    " << res[i][6] << endl;
    }*/


    fstream file("file.txt");
    for (int i=0; i < cnt; i++) {
        file <<  res[i][1] << " " << res[i][2] << " " << res[i][3]  << endl;
    }

    file.close();


    fstream file_station("stations.txt");
    for (int i=0; i < cnt; i++) {
        double** stations = create_observatories(res[i][0]);
        double distance = pow(res[i][1],2) + pow(res[i][2],2) + pow(res[i][3],2);
        double max_distance = sqrt(distance - pow(R_CONST,2));
        for (int j = 0; j < 8; j++) {
            double r = sqrt(pow((stations[j][0] - res[i][1]), 2) + pow((stations[j][1] - res[i][2]), 2) +
                              pow((stations[j][2] - res[i][3]), 2));
            if (r <= max_distance){
               file_station << res[i][0] - JD_start << " " << j << " "<< r << endl;
               //cout << res[i][0] - JD_start << " " << j << " "<< r << endl;

            }
        }
    }
    file_station.close();


    for(int i=0; i < cnt; i++){
        delete[] res[i];
    }
    delete[] res;

    delete[] vec;

}
#include "main.h"

#define JD 2451545.0
#define STEP 30.0
#define GM 398600.4415 // км^3 / с^2
#define START_POINT 6878.0 //км
#define START_POINT_VAR 6878.0 //км
#define GM_VAR 398601.45 // км^3/с^2
#define J2_VAR 1.75553e10


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


    double *states = new double[54];
    states[0] = START_POINT_VAR, states[1] = 0, states[2] = 0;
    states[3] = 0, states[4] = sqrt(GM_VAR / START_POINT_VAR), states[5] = 0;
    for(int i=6; i < 54; i++){
        states[i] = 0;
    }
    states[6] = 1; states[13] = 1; states[20] = 1; states[27] = 1; states[34] = 1; states[41] = 1;


    iauC2t06a(JD_start + (37.0 + 32.184) / 86400.0, 0, JD_start, 0, 0, 0, rotateMatrix);
    Transposition(rotateMatrix);
    changeCoords(rotateMatrix, vec, 0);
    changeCoords(rotateMatrix, vec, 3);


    for (int i=0; i < 2; i++){ ///?
        changeCoords(rotateMatrix, states, 3*i);  //НСК
    }

    double *b  =new double[8];
    for (int i=0; i < 6; i++){
        b[i] = states[i];
    }
    b[6] = GM_VAR;
    b[7] = J2_VAR;


    double** new_res = integrate_for_inverse(JD, STEP, 54, states, J2_VAR, GM_VAR);
    double** res = integrate(JD, STEP, 6, vec);

    /*
    for (int i=0; i < cnt; i++){
        for (int j=0; j < 55; j++){
            cout << new_res[i][j] << " ";
        }
        cout << endl << endl;
    }
    cout << endl << endl; */

    vector<vector<double>> A;

    vector<double> r_b;


    for (int i=0; i < 2 * cnt; i++) {
        double **stations = create_observatories(res[i][0]);
        double distance = pow(res[i][1], 2) + pow(res[i][2], 2) + pow(res[i][3], 2);
        double max_distance = sqrt(distance - pow(R_CONST, 2));

        for(int station_number = 0; station_number < 8; station_number++) {
            double r_var = sqrt(pow((stations[station_number][0] - new_res[i][1]), 2) +
                                pow((stations[station_number][1] - new_res[i][2]), 2) +
                                pow((stations[station_number][2] - new_res[i][3]), 2));
            double r_original = sqrt(
                    pow((stations[station_number][0] - res[i][1]), 2) +
                    pow((stations[station_number][1] - res[i][2]), 2) +
                    pow((stations[station_number][2] - res[i][3]), 2));

            if (r_original <= max_distance) {
                double *dg_dX = new double[6];
                for (int t = 3; t < 6; t++) {
                    dg_dX[t] = 0;
                }
                dg_dX[0] = dg_dx(new_res[i], stations[station_number]); ///
                dg_dX[1] = dg_dy(new_res[i], stations[station_number]);
                dg_dX[2] = dg_dz(new_res[i], stations[station_number]);

                vector<double> current_r;

                for (int k = 0; k < 8; k++) {
                    double result = 0;
                    for (int t = 0; t < 6; t++) {
                        result += dg_dX[t] * new_res[i][7 + 6 * k + t];
                    }
                    current_r.push_back(-result);
                }
                A.push_back(current_r);
                r_b.push_back(r_var - r_original);
            }

        }
    }


    double** AtA = multiplication_AtA(A);
    double* Atr = multiplication_Atr(A, r_b);


    double* x = Cholesky_decomposition(AtA, 8, Atr);

    double *b_new = new double[8];


    for(int i=0; i < 8; i++){
        b_new[i] = b[i] - x[i];
    }


    for(int i=0; i < 8; i++){
        cout << b[i] << "  ";
    }
    cout << endl;

    for(int i=0; i < 8; i++){
        cout << b_new[i] << "  ";
    }
    cout << endl;



    fstream file("file.txt");
    for (int i=0; i < cnt; i++) {
        file <<  res[i][1] << " " << res[i][2] << " " << res[i][3]  << endl;
    }

    file.close();





    for(int i=0; i < 2* cnt; i++){
        delete[] res[i];
        delete[] new_res[i];
    }
    delete[] res;
    delete[] new_res;

    delete[] vec;
    delete[] states;


}
#include "main.h"

#define JD 2451545.0
#define STEP 60.0
#define GM 398600.4415 // км^3 / с^2
#define START_POINT 6878.0 //км
#define J2 1.75553e10



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

    double noise[8] = {0.0005,-0.0005,0.0002,0.0000005,-0.0000005,0.00000003,2, 10};

    int cnt = GENERAL_TIME / STEP;
    double JD_start = JD;
    double *vec = new double[6];
    double rotateMatrix[3][3];
    vec[0] = START_POINT, vec[1] = 0, vec[2] = 0;
    vec[3] = 0, vec[4] = sqrt(GM / START_POINT), vec[5] = 0;


    double *states = new double[54];
    states[0] = START_POINT, states[1] = 0, states[2] = 0;
    states[3] = 0, states[4] = sqrt(GM / START_POINT), states[5] = 0;
    for(int i=6; i < 54; i++){
        states[i] = 0;
    }
    states[6] = 1; states[13] = 1; states[20] = 1; states[27] = 1; states[34] = 1; states[41] = 1;

    double *init = new double[8];
    for (int i=0; i < 6; i++){
        init[i] = states[i];
    }
    init[6] = GM;
    init[7] = J2;

    for (int i=0; i < 6; i++){
        states[i] += noise[i];
    }

    iauC2t06a(JD_start + (37.0 + 32.184) / 86400.0, 0, JD_start, 0, 0, 0, rotateMatrix);
    Transposition(rotateMatrix);
    changeCoords(rotateMatrix, vec, 0);
    changeCoords(rotateMatrix, vec, 3);
    changeCoords(rotateMatrix, init, 0);
    changeCoords(rotateMatrix, init, 3);

    for (int i=0; i < 8; i++){
        cout << init[i] << "  ";
    }
    cout << endl << endl;


    for (int i=0; i < 18; i++){ ///?
        changeCoords(rotateMatrix, states, 3*i);  //НСК
    }


    double *b = new double[8];
    for (int i=0; i < 6; i++){
        b[i] = states[i];
    }

    b[6] = GM + noise[6];
    b[7] = J2 + noise[7];


    double** res = integrate(JD, STEP, 6, vec);


    for(int iter = 0; iter < 4; iter++) {
        double **new_res = integrate_for_inverse(JD, STEP, 54, states, b[7], b[6]);
        vector<vector<double>> A;
        vector<double> r_b;

        for (int i = 0; i <  cnt; i++) {
            double **stations = create_observatories(res[i][0]);
            double distance = pow(res[i][1], 2) + pow(res[i][2], 2) + pow(res[i][3], 2);
            double max_distance = sqrt(distance - pow(R_CONST, 2));

            for (int station_number = 0; station_number < 8; station_number++) {
                double r_var = sqrt(pow((stations[station_number][0] - new_res[i][1]), 2) +
                                    pow((stations[station_number][1] - new_res[i][2]), 2) +
                                    pow((stations[station_number][2] - new_res[i][3]), 2));
                double r_original = sqrt(
                        pow((stations[station_number][0] - res[i][1]), 2) +
                        pow((stations[station_number][1] - res[i][2]), 2) +
                        pow((stations[station_number][2] - res[i][3]), 2));

                double r_original_var = r_original * 1.0;

                if (r_original <= max_distance) {
                    double *dg_dX = new double[6];
                    for (int t = 3; t < 6; t++) {
                        dg_dX[t] = 0;
                    }
                    dg_dX[0] = dg_dx(new_res[i], stations[station_number]);
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
                    /*
                    for (int r=0; r < 8; r++){
                        cout << current_r[r] << " ";
                    }
                    cout << endl; */

                    A.push_back(current_r);
                    r_b.push_back(r_var - r_original_var);
                }

            }
        }

        double var = 1.0/(pow(0.00001, 2));

        double **AtA = multiplication_AtA(A);

        for (int i=0; i < 8; i++){
            for (int j=0; j < 8; j++){
                AtA[i][j] *= var;
            }
        }

        double *Atr = multiplication_Atr(A, r_b);
        for (int i=0; i < 8; i++){
            Atr[i] *= var;
        }


        double *x = Cholesky_decomposition(AtA, 8, Atr);

        double *b_new = new double[8];


        for (int i = 0; i < 8; i++) {
            b_new[i] = b[i] - x[i];
        }

        for (int i = 0; i < 8; i++) {
            cout << b[i] << " | " << b_new[i] << "  ";
        }
        cout << endl;

        for (int i = 0; i < 8; i++) {
            b[i] = b_new[i];
        }

        for (int i = 0; i < 6; i++){
            states[i] = b[i];
        }

        for(int i = 6; i < 54; i++){
            states[i] = 0;
        }

        states[6] = 1; states[13] = 1; states[20] = 1; states[27] = 1; states[34] = 1; states[41] = 1;

    }



    fstream file("file.txt");
    for (int i=0; i < cnt; i++) {
        file <<  res[i][1] << " " << res[i][2] << " " << res[i][3]  << endl;
    }

    file.close();


    for(int i=0; i < cnt; i++){
        delete[] res[i];

    }
    delete[] res;


    delete[] vec;
    delete[] states;


}

#include "main.h"

#define JD 2451545.0
#define STEP 60.0
#define GM 398600.4415 // км^3 / с^2
#define START_POINT 6878.0 //км
#define J2 1.75553e10
#define STATIONS_NUMBER 4



double** create_observatories(double JD_start){

    Observatories initial_coord[18] = {
            Observatories(281.5075, 1.00045, -0.00405),
            Observatories(260.8053, 0.94388, +0.33026),
            Observatories(107.6160, 0.99316, -0.11808),
            Observatories(281.65, 0.999, +0.000),
            Observatories(321.3126, 0.98840, -0.15179),
            Observatories(299.99039, 0.998647, -0.051941),
            Observatories(101.43942, 0.998617, +0.052565),
            Observatories(282.70896, 1.000183, +0.021030),
            Observatories(98.48553, 0.948521, +0.316891),
            Observatories(99.78111, 0.994005, +0.109127),
            Observatories(101.27869, 0.975556, +0.218996),
            Observatories(114.08987, 0.928304, -0.370597),
            Observatories(203.74299, 0.936235, +0.351547),
            Observatories(210.39020, 0.953686, -0.299837),
            Observatories(204.52398, 0.941706, +0.337237),
            Observatories(286.36626, 0.995574, +0.097155),
            Observatories(289.32156, 0.957949, -0.287797),
            Observatories( 289.32156, 0.957949, -0.287797),

    };

    double** station = new double *[STATIONS_NUMBER];
    for (int i=0; i < STATIONS_NUMBER; i++){
        station[i] = initial_coord[i].get_coords();
    }

    double rotateMatrix[3][3];
    iauC2t06a(JD_start + (37.0 + 32.184) / 86400.0, 0, JD_start, 0, 0, 0, rotateMatrix);
    Transposition(rotateMatrix);
    for (int i=0; i < STATIONS_NUMBER; i++){
        changeCoords(rotateMatrix, station[i], 0);
    }

    return station;

}

double** create_model_observatories(double JD_start){
    double** station = new double *[STATIONS_NUMBER];
    for (int i=0; i < STATIONS_NUMBER; i++){
        station[i] = new double[3];
    }
    station[0][0] = 6378.1363 ; station[0][1] = 0; station[0][2] = 0;
    station[1][0] = -6378.1363; station[1][1] = 0; station[1][2] = 0;
    station[2][0] = 0; station[2][1] = 6378.1363; station[2][2] = 0;
    station[3][0] = 0; station[3][1] = -6378.1363; station[3][2] = 0;


    double rotateMatrix[3][3];
    iauC2t06a(JD_start + (37.0 + 32.184) / 86400.0, 0, JD_start, 0, 0, 0, rotateMatrix);
    Transposition(rotateMatrix);
    for (int i=0; i < STATIONS_NUMBER; i++){
        changeCoords(rotateMatrix, station[i], 0);
    }

    return station;
}

int main(){

    double noise[8] = {0.0005,-0.0005,0.0000002,0.0000005,-0.0000005,0.0000002,2, -10};
    //double noise[8] = {0,0,0,0,0,0,0, 0};
    //double noise[8] = {0.0005,-0.0005,0,0.0000005,-0.0000005,0,2, -10};

    int cnt = GENERAL_TIME  / STEP;
    double JD_start = JD;
    double *vec = new double[6];

    vec[0] = 1248.77, vec[1] = -6763.69, vec[2] = -0.155766;
    vec[3] = 7.48616, vec[4] = 1.38216, vec[5] = 0.00024043;

    //vec[0] = 1472.62, vec[1] = -6718.5, vec[2] = -0.148523;
    //vec[3] = 7.43608, vec[4] = 1.63027, vec[5] = 0.000242242;

    //vec[0] = START_POINT, vec[1] = 0, vec[2] = 0;
    //vec[3] = 0, vec[4] = sqrt(GM / START_POINT), vec[5] = 0;

    for (int i=0; i < 6; i++){
        cout << vec[i] << "  ";
    }
    cout << GM << "  " << J2;
    cout << endl << endl;

    double *states = new double[54];
    for(int i=0; i < 6; i++){
        states[i] = vec[i];
    }
    for(int i=6; i < 54; i++){
        states[i] = 0;
    }
    states[6] = 1; states[13] = 1; states[20] = 1; states[27] = 1; states[34] = 1; states[41] = 1;

    for (int i=0; i < 6; i++){
        states[i] += noise[i];
    }


    double *b = new double[8];
    for (int i=0; i < 6; i++){
        b[i] = states[i];
    }

    b[6] = GM + noise[6];
    b[7] = J2 + noise[7];


    double** res = integrate(JD, STEP, 6, vec);
    /*
    for (int i=0; i < cnt; i++ ){
        for (int j=0; j < 7; j++){
            cout << res[i][j] << " ";
        }
        cout << endl;
    }*/

    double **new_res;


    /*
    random_device randomDevice;
    mt19937 generation(randomDevice());
    double var = 1.0/(pow(0.00001, 2));
    uniform_real_distribution<double> distribution(-0.5, 0.5);
    double W[8];
    for (int i=0; i < 8; i++){
        W[i] = 1.0 / pow( 1.0/ 10000, 2)   ;
    }

    for (int i=0; i < 8; i++){
        cout << W[i] << " ";
    }
    cout << endl;

    */

    for(int iter = 0; iter < 4; iter++) {

        new_res = integrate_for_inverse(JD, STEP, 54, states, b[7], b[6]);

        vector<vector<double>> A;
        vector<double> r_b;
        random_device rd;
        mt19937 gen(rd());

        for (int i = 0; i <  cnt; i++) {

            double **stations = create_model_observatories(res[i][0]);

            double distance = pow(res[i][1], 2) + pow(res[i][2], 2) + pow(res[i][3], 2);
            double max_distance = sqrt(distance - pow(R_CONST, 2));

            for (int station_number = 0 ; station_number < STATIONS_NUMBER; station_number++) {

                double r_var = sqrt(pow((stations[station_number][0] - new_res[i][1]), 2) +
                                    pow((stations[station_number][1] - new_res[i][2]), 2) +
                                    pow((stations[station_number][2] - new_res[i][3]), 2));
                double r_original = sqrt(
                        pow((stations[station_number][0] - res[i][1]), 2) +
                        pow((stations[station_number][1] - res[i][2]), 2) +
                        pow((stations[station_number][2] - res[i][3]), 2));

                //cout << distance << " " << r_original << endl;

                uniform_real_distribution<double> dist(-r_original * 0.01, r_original * 0.01);
                double r_original_var = r_original ;//+ dist(gen);

                if (r_original <= max_distance) {
                    //cout << station_number << " ";

                    double *dg_dX = new double[6];
                    for (int w = 3; w < 6; w++) {
                        dg_dX[w] = 0;
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

                    A.push_back(current_r);
                    r_b.push_back(r_var - r_original_var);
                }

            }
        }


        double **AtA = multiplication_AtA(A);
        double *Atr = multiplication_Atr(A, r_b);

        /*
        double W[8];
        for (int i=0; i < 8; i++){
            W[i] = 1.0 / pow( 0.00001, 2)   ;
        }

        for (int i=0; i < 8; i++){
            for (int j=0; j < 8; j++){
                AtA[i][j] *= W[i];
            }
        }

        for (int i=0; i < 8; i++){
            Atr[i] *= W[i];
        } */

        double *x = Cholesky_decomposition(AtA, 8, Atr);
        /*cout << "X: ";
        for (int t=0; t < 8; t++){
            cout << x[t] << " ";
        }
        cout << endl; */

        double *b_new = new double[8];


        for (int i = 0; i < 8; i++) {
            b_new[i] = -x[i] + b[i];
        }
        cout << endl;
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
        file <<  res[i][1] << " " << res[i][2] << " " << res[i][3] << endl;
    }

    file.close();


    for(int i=0; i < cnt; i++){
        delete[] res[i];

    }
    delete[] res;


    delete[] vec;
    delete[] states;


}
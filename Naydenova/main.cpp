#include "main.h"

#define JD 2451545.0
#define STEP 60.0
#define GM 398600.4415 // км^3 / с^2
#define START_POINT 6878.0 //км
#define J2 1.75553e10
#define STATIONS_NUMBER 4



double** create_observatories(double JD_start){

    Observatories initial_coord[80] = {
            Observatories(55.5061, 0.93464, -0.35447), // St. Clotilde, Reunion
            Observatories(55.4100, 0.93288, -0.35941), // Observatoire des Makes, Saint-Louis
            Observatories(55.2586, 0.93394, -0.35634), // St. Paul, Reunion
            Observatories(78.4541, 0.95444, +0.29768), // Hyderabad
            Observatories(78.7283, 0.95618, +0.29216), // Japal-Rangapur
            Observatories(78.8263, 0.97627, +0.21634), // Vainu Bappu Observatory, Kavalur
            Observatories(80.2464, 0.97427, +0.22465), // Madras
            Observatories(293.24692, 0.949577, +0.312734), // Arecibo
            Observatories(204.52344, 0.941701, +0.337237), // New Horizons KBO Search-Subaru
            Observatories(204.53044, 0.941705, +0.337234), // New Horizons KBO Search-CFHT
            Observatories(107.6160, 0.99316, -0.11808), // Bosscha Observatory, Lembang
            Observatories(288.88, 0.990, +0.150), // University of the Andes station
            Observatories(289.1296, 0.98890, +0.15185), // OAN de Llano del Hato, Merida
            Observatories(290.6769, 0.98477, +0.17381), // Observatorio Taya Beixo, Barquisimeto
            Observatories(112.334, 0.9574, +0.2877), // Tsingtao field station, Xisha Islands
            Observatories(203.7424, 0.93623, +0.35156), // Haleakala-NEAT/GEODSS
            Observatories(204.5278, 0.94171, +0.33725), // Maunakea
            Observatories(203.7420, 0.93623, +0.35156), // Haleakala-AMOS
            Observatories(203.5683, 0.93557, +0.35201), // Kihei-AMOS Remote Maui Experimental Site
            Observatories(260.8053, 0.94388, +0.33026), // National Observatory, Tacubaya
            Observatories(263.2300, 0.95591, +0.29359), // Oaxaca
            Observatories(281.5075, 1.00045, -0.00405), // Quito
            Observatories(281.65, 0.999, +0.000), // Quito, comet astrograph station
            Observatories(288.4511, 0.96006, -0.28021), // Harvard Observatory, Arequipa
            Observatories(295.37581, 0.930491, -0.365872), // Tarija
            Observatories(316.3097, 0.94132, -0.33707), // Wykrota Observatory-CEAMIG
            Observatories(210.41224, 0.953752, -0.299638), // S. S. Observatory, Pamatai
            Observatories(210.3842, 0.95330, -0.30100), // Puna'auia
            Observatories(120.6982, 0.92730, +0.37308), // Kenting Observatory, Checheng
            Observatories(120.7839, 0.92796, +0.37148), // Kenting Observatory, Hengchun
            Observatories(145.7403, 0.96545, +0.25977), // Pacific Sky Observatory, Saipan
            Observatories(145.69721, 0.957625, -0.287092), // Earl Hill Observatory, Trinity Beach
            Observatories(203.74409, 0.936241, +0.351543), // Pan-STARRS 1, Haleakala
            Observatories(203.74409, 0.936239, +0.351545), // Pan-STARRS 2, Haleakala
            Observatories(201.94100, 0.932037, +0.361160), // Ironwood Remote Observatory, Hawaii
            Observatories(201.95283, 0.929942, +0.366558), // Ironwood Observatory, Hawaii
            Observatories(203.74250, 0.936239, +0.351538), // Haleakala-Faulkes Telescope North
            Observatories(210.3842, 0.95330, -0.30100), // Hibiscus Observatory, Punaauia
            Observatories(210.3842, 0.95330, -0.30100), // Tiki Observatory, Punaauia
            Observatories(210.38381, 0.953304, -0.301004), // Moana Observatory, Punaauia
            Observatories(282.70896, 1.000183, +0.021030), // University of Narino Observatory, Pasto
            Observatories(312.0000, 0.96235, -0.27153), // Taurus Australis Observatory, Brasilia
            Observatories(312.4981, 0.96931, -0.24573), // Pousada dos Anoes Observatory
            Observatories(316.0025, 0.94119, -0.33714), // CEAMIG-REA Observatory, Belo Horizonte
            Observatories(39.25827, 0.930706, +0.364567), // Jeddah
            Observatories(98.48553, 0.948521, +0.316891), // TRT-NEO, Chiangmai
            Observatories(99.78111, 0.994005, +0.109127), // Observatori Negara, Langkawi
            Observatories(101.43942, 0.998617, +0.052565), // Hin Hua Observatory, Klang
            Observatories(101.27869, 0.975556, +0.218996), // Akin Observatory, Rayong
            Observatories(114.08987, 0.928304, -0.370597), // Space Surveillance Telescope, HEH Station
            Observatories(203.74247, 0.936240, +0.351538), // Haleakala-LCO Clamshell #3
            Observatories(203.74249, 0.936241, +0.351538), // Haleakala-LCO OGG B #2
            Observatories(203.74299, 0.936235, +0.351547), // ATLAS-HKO, Haleakala
            Observatories(204.42387, 0.943290, +0.332467), // ATLAS-MLO Auxiliary Camera, Mauna Loa
            Observatories(204.42395, 0.943290, +0.332467), // ATLAS-MLO, Mauna Loa
            Observatories(204.52398, 0.941706, +0.337237), // Subaru Telescope, Maunakea
            Observatories(204.52241, 0.941706, +0.337212), // Submillimeter Array, Maunakea (SMA)
            Observatories(204.53036, 0.941731, +0.337198), // United Kingdom Infrared Telescope, Maunakea
            Observatories(204.53057, 0.941729, +0.337199), // University of Hawaii 88-inch telescope, Maunakea
            Observatories(204.52771, 0.941691, +0.337263), // NASA Infrared Telescope Facility, Maunakea
            Observatories(204.53113, 0.941714, +0.337236), // Canada-France-Hawaii Telescope, Maunakea
            Observatories(204.53094, 0.941727, +0.337214), // Gemini North Observatory, Maunakea
            Observatories(204.52570, 0.941703, +0.337250), // W. M. Keck Observatory, Keck 1, Maunakea
            Observatories(204.52580, 0.941700, +0.337256), // W. M. Keck Observatory, Keck 2, Maunakea
            Observatories(210.39020, 0.953686, -0.299837), // Astronomical Society of Tahiti
            Observatories(259.69219, 0.936737, +0.349756), // Observatoire LAURIER, El Marques
            Observatories(286.36626, 0.995574, +0.097155), // AstroExplor Observatory, Tinjaca
            Observatories(284.30958, 0.996767, +0.082976), // Observatorio Astronomico UTP, Pereira
            Observatories(289.32156, 0.957949, -0.287797), // Observatorio Astronomico de Moquegua, Carumas
            Observatories(299.99039, 0.998647, -0.051941), // OARU, Manaus
            Observatories(308.43322, 0.931623, -0.362375), // Observatorio OATU, Tupi Paulista
            Observatories(309.1506, 0.92987, -0.36684), // Observatorio Campo dos Amarais
            Observatories(312.0889, 0.96218, -0.27210), // Dogsheaven Observatory, Brasilia
            Observatories(312.21786, 0.962401, -0.271311), // Rocca Observatory, Brasilia
            Observatories(312.13208, 0.963322, -0.268189), // Carina Observatory, Brasilia
            Observatories(315.21504, 0.935906, -0.351562), // SONEAR Observatory, Oliveira
            Observatories(316.01580, 0.940890, -0.337985), // SONEAR 2 Observatory, Belo Horizonte
            Observatories(318.68794, 0.929268, -0.368170), // ROCG, Campos dos Goytacazes
            Observatories(321.3126, 0.98840, -0.15179), // OASI, Nova Itacuruba
            Observatories(324.03889, 0.989706, -0.143217), // Discovery Observatory, Caruaru
    };

    double** station = new double *[STATIONS_NUMBER];
    for (int i=0; i < STATIONS_NUMBER; i++){
        station[i] = initial_coord[i].get_coords();
    }

    double rotateMatrix[3][3];
    rotateMatrix[0][0] = 1;
    rotateMatrix[1][1] = 1;
    rotateMatrix[2][2] = 1;
    //iauC2t06a(JD_start + (37.0 + 32.184) / 86400.0, 0, JD_start, 0, 0, 0, rotateMatrix);
    Transposition(rotateMatrix);
    for (int i=0; i < STATIONS_NUMBER; i++){
        changeCoords(rotateMatrix, station[i], 0);
    }

    return station;

}

int main(){

    double noise[8] = {0.0005,-0.0005,0.0000002,0.0000005,-0.0000005,0.0000002,2, -10};
    //double noise[8] = {0.0005,-0.0005,0,0.0000005,-0.0000005,0,2, -10};

    int cnt = GENERAL_TIME / STEP;
    double JD_start = JD;
    double *vec = new double[6];

    //vec[0] = 1248.77, vec[1] = -6763.69, vec[2] = -0.155766;
    //vec[3] = 7.48616, vec[4] = 1.38216, vec[5] = 0.00024043;
    vec[0] = START_POINT, vec[1] = 0, vec[2] = 0;
    vec[3] = 0, vec[4] = sqrt(GM / START_POINT), vec[5] = 0;

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

    double ** stations = new double* [STATIONS_NUMBER];
    for (int i=0; i < STATIONS_NUMBER; i++){
        stations[i] = new double [3];
    }
    stations[0][0] = 6378.1363 ; stations[0][1] = 0; stations[0][2] = 0;
    stations[1][0] = -6378.1363; stations[1][1] = 0; stations[1][2] = 0;
    stations[2][0] = 0; stations[2][1] = 6378.1363; stations[2][2] = 0;
    stations[3][0] = 0; stations[3][1] = -6378.1363; stations[3][2] = 0;


    for(int iter = 0; iter < 4; iter++) {

        new_res = integrate_for_inverse(JD, STEP, 54, states, b[7], b[6]);
        /*
        for (int i=0; i < cnt; i++){
            for (int j=0; j < 55; j++){
                cout << new_res[i][j] << " ";
            }
            cout << endl << endl;
        } */

        vector<vector<double>> A;
        vector<double> r_b;
        random_device rd;
        mt19937 gen(rd());

        for (int i = 0; i < cnt; i++) {

            //double **stations = create_observatories(res[i][0]);

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

                uniform_real_distribution<double> dist(-r_original * 0.02, r_original * 0.02);
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

        //double var = 1.0/(pow(0.00001, 2));

        double **AtA = multiplication_AtA(A);

        /*
        for (int i=0; i < 8; i++){
            for (int j=0; j < 8; j++){
                AtA[i][j] *= var;
            }
        } */

        double *Atr = multiplication_Atr(A, r_b);
        /*
        for (int i=0; i < 8; i++){
            Atr[i] *= var;
        } */


        double *x = Cholesky_decomposition(AtA, 8, Atr);
        /*cout << "X: ";
        for (int t=0; t < 8; t++){
            cout << x[t] << " ";
        }
        cout << endl; */

        double *b_new = new double[8];


        for (int i = 0; i < 8; i++) {
            b_new[i] = b[i] - x[i];
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
#include "AdvancedFunction.h"

double** create_matrix_df_dx(double* x, double mu, double J){
    double **df_dx = new double*[6];
    for(int i=0; i < 6;i++){
        df_dx[i] = new double[6];
        for (int j = 0; j < 6; j++){
            df_dx[i][j] = 0;
        }
    }

    df_dx[0][3] = 1;
    df_dx[1][4] = 1;
    df_dx[2][5] = 1;

    df_dx[3][0] = -dux_dx(x, mu, J);
    df_dx[3][1] = -dux_dy(x, mu, J);
    df_dx[3][2] = -dux_dz(x, mu, J);
    df_dx[4][0] = -duy_dx(x, mu, J);
    df_dx[4][1] = -duy_dy(x, mu, J);
    df_dx[4][2] = -duy_dz(x, mu, J);
    df_dx[5][0] = -duz_dx(x, mu, J);
    df_dx[5][1] = -duz_dy(x, mu, J);
    df_dx[5][2] = -duz_dz(x, mu, J);

    return df_dx;
}

void function(double* x, double* vec, double JD, double J, double mu){
    vec[0] = x[3];
    vec[1] = x[4];
    vec[2] = x[5];
    vec[3] = 0;
    vec[4] = 0;
    vec[5] = 0;


    double rotateMatrix[3][3];

    iauC2t06a(JD + (37.0 + 32.184) / 86400.0, 0, JD, 0, 0, 0, rotateMatrix);

    changeCoords(rotateMatrix, x, 0);

    double *grad = new double[3];
    grad[0] = -dx(x);
    grad[1] = -dy(x);
    grad[2] = -dz(x);

    double *states = new double[48];

    double **df_dx = create_matrix_df_dx(x, mu, J);


   for (int i=0; i < 6; i++){
       for (int j = 0; j < 8; j++){
           double res = 0;
           for (int k = 0; k < 6; k++){
               res += df_dx[i][k] * x[6 + 6*j + k];
           }
           states[6*j + i] = res;
       }
   }

   states[39] += -dux_dmu(x);
   states[40] += -duy_dmu(x);
   states[41] += -duz_dmu(x);
   states[45] += -dux_dJ(x);
   states[46] += -duy_dJ(x);
   states[47] += -duz_dJ(x);

    Transposition(rotateMatrix);

    changeCoords(rotateMatrix, grad, 0);

    for (int i=0; i < 16; i++){  ///?
        changeCoords(rotateMatrix, states, 3*i);
    }

    for (int i = 0; i < 3; i++) {
        vec[i + 3] = grad[i];
    }

    for(int i= 0; i < 48; i++){
        vec[i + 6] = states[i];
    }
}
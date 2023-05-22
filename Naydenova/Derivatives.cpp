#include "Derivatives.h"


double dx(double* vec) {

    double sum_sqrs = pow(vec[0],2) + pow(vec[1],2) + pow(vec[2],2);
    double res = MU_CONST * (vec[0] / pow(sum_sqrs, 1.5)) +
            J2 * 1.5 * (( (-5*vec[0]*pow(vec[2],2))/pow(sum_sqrs, 3.5))
            - (1.0/3.0)* ( -3*vec[0] / pow(sum_sqrs, 2.5) ));
    return res;
}

double dy(double* vec) {

    double sum_sqrs = pow(vec[0],2) + pow(vec[1],2) + pow(vec[2],2);
    double res = MU_CONST * (vec[1] / pow(sum_sqrs, 1.5)) +
                 J2 * 1.5 * (( (-5*vec[1]*pow(vec[2],2))/pow(sum_sqrs, 3.5))
                 - (1.0/3.0)* ( -3*vec[1] / pow(sum_sqrs, 2.5) ));

    return res;
}

double dz(double* vec) {
    double sum_sqrs = pow(vec[0],2) + pow(vec[1],2) + pow(vec[2],2);
    double res = MU_CONST * (vec[2] / pow(sum_sqrs, 1.5)) +
                 J2 * 1.5 * (((vec[2]* (2 * pow(vec[0],2) + 2 * pow(vec[1],2) - 3 * pow(vec[2],2)))/pow(sum_sqrs, 3.5))
                 - (1.0/3.0)* ( -3*vec[2] / pow(sum_sqrs, 2.5) ));

    return res;
}

double dg_dx(double* vec, double* coord){

    double res = (vec[1] - coord[0])/( sqrt(pow(vec[1] - coord[0], 2) + pow(vec[2] - coord[1], 2) + pow(vec[3] - coord[2], 2) ));
    return res;
}

double dg_dy(double* vec, double* coord){
    double res = (vec[2] - coord[1])/( sqrt(pow(vec[1] - coord[0], 2) + pow(vec[2] - coord[1], 2) + pow(vec[3] - coord[2], 2) ));
    return res;
}

double dg_dz(double* vec, double* coord){
    double res = (vec[3] - coord[2])/( sqrt(pow(vec[1] - coord[0], 2) + pow(vec[2] - coord[1], 2) + pow(vec[3] - coord[2], 2) ));
    return res;
}
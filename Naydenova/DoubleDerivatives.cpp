#include "DoubleDerivatives.h"


double dux_dx(double* vec, double mu, double J){
    double sum_sqrs = pow(vec[0],2) + pow(vec[1],2) + pow(vec[2],2);
    double res = -mu * ((2* pow(vec[0],2) - pow(vec[1],2) - pow(vec[2], 2))/pow(sum_sqrs,2.5)) +
                 J * 1.5 * (((35 * pow(vec[0],2) *pow(vec[2],2) * pow(sum_sqrs, 0.5) *
                 (pow(vec[0],4) + pow(vec[1],4) + pow(vec[2],4) + 2*pow(vec[0],2)*pow(vec[1],2) + 2*pow(vec[0],2)*pow(vec[2],2)
                 + 2*pow(vec[1],2)*pow(vec[2],2)) - 5*pow(vec[2],2)* pow(sum_sqrs,3.5)) / pow(sum_sqrs,7))
                 - (1.0/3.0)* ((15*pow(vec[0],2) * pow(sum_sqrs,0.5) - 3* pow(sum_sqrs, 1.5)) / pow(sum_sqrs,4)));
    return res;
}

double dux_dy(double* vec, double mu, double J){
    double sum_sqrs = pow(vec[0],2) + pow(vec[1],2) + pow(vec[2],2);
    double res = -mu * ((3*vec[0] * vec[1])/pow(sum_sqrs,2.5)) +
                 J * 1.5 * (((35*vec[0]* vec[1] * pow(vec[2],2))/ (pow(sum_sqrs,4.5)))
                 - (1.0/3.0)* ((15*vec[0]*vec[1] ) / pow(sum_sqrs,3.5)));
    return res;
}

double dux_dz(double* vec, double mu, double J){
    double sum_sqrs = pow(vec[0],2) + pow(vec[1],2) + pow(vec[2],2);
    double res = -mu * ((3*vec[0] * vec[2])/pow(sum_sqrs,2.5))  +
                 J * 1.5 * (((35 * vec[0] *pow(vec[2],3) * pow(sum_sqrs,0.5) *
                 (pow(vec[0],4) + pow(vec[1],4) + pow(vec[2],4) + 2*pow(vec[0],2)*pow(vec[1],2) + 2*pow(vec[0],2)*pow(vec[2],2)
                 + 2*pow(vec[1],2)*pow(vec[2],2)) - 10 * vec[0] * vec[2] * pow(sum_sqrs,3.5)) / pow(sum_sqrs,7))
                 - (1.0/3.0)* ((15*vec[0]*vec[2] ) / pow(sum_sqrs,3.5)));
    return res;
}

double duy_dx(double* vec, double mu, double J){
    double sum_sqrs = pow(vec[0],2) + pow(vec[1],2) + pow(vec[2],2);
    double res = -mu * ((3*vec[0] * vec[1])/pow(sum_sqrs,2.5)) +
                 J * 1.5 * (((35*vec[0]* vec[1] * pow(vec[2],2))/ (pow(sum_sqrs,4.5)))
                 - (1.0/3.0)* ((15*vec[0]*vec[1] ) / pow(sum_sqrs,3.5)));
    return res;
}


double duy_dy(double* vec, double mu, double J){
    double sum_sqrs = pow(vec[0],2) + pow(vec[1],2) + pow(vec[2],2);
    double res = -mu * ((2* pow(vec[1],2) - pow(vec[0],2) - pow(vec[2], 2))/pow(sum_sqrs,2.5)) +
                 J * 1.5 * (((35 * pow(vec[1],2) *pow(vec[2],2) * pow(sum_sqrs, 0.5) *
                 (pow(vec[0],4) + pow(vec[1],4) + pow(vec[2],4) + 2*pow(vec[0],2)*pow(vec[1],2) + 2*pow(vec[0],2)*pow(vec[2],2)
                 + 2*pow(vec[1],2)*pow(vec[2],2)) - 5*pow(vec[2],2)* pow(sum_sqrs,3.5)) / pow(sum_sqrs,7))
                 - (1.0/3.0)* ((15*pow(vec[1],2) * pow(sum_sqrs,0.5) - 3* pow(sum_sqrs, 1.5)) / pow(sum_sqrs,4)));
    return res;
}


double duy_dz(double* vec, double mu, double J){
    double sum_sqrs = pow(vec[0],2) + pow(vec[1],2) + pow(vec[2],2);
    double res = -mu * ((3*vec[1] * vec[2])/pow(sum_sqrs,2.5))  +
                 J * 1.5 * (((35 * vec[1] *pow(vec[2],3) * pow(sum_sqrs,0.5) *
                 (pow(vec[0],4) + pow(vec[1],4) + pow(vec[2],4) + 2*pow(vec[0],2)*pow(vec[1],2) + 2*pow(vec[0],2)*pow(vec[2],2)
                 + 2*pow(vec[1],2)*pow(vec[2],2)) - 10 * vec[1] * vec[2] * pow(sum_sqrs,3.5)) / pow(sum_sqrs,7))
                 - (1.0/3.0)* ((15*vec[1]*vec[2] ) / pow(sum_sqrs,3.5)));
    return res;
}


double duz_dx(double* vec, double mu, double J){
    double sum_sqrs = pow(vec[0],2) + pow(vec[1],2) + pow(vec[2],2);
    double res = -mu * ((3*vec[0] * vec[2])/pow(sum_sqrs,2.5)) +
                 J * 1.5 * (((-8 * pow(vec[0], 3) * vec[2] - 8 * vec[0] * pow(vec[1], 2) * vec[2] + 22 * vec[0] * pow(vec[2], 3)) / pow(sum_sqrs, 4))
                            - (1.0/3.0)* ((15*vec[0]*vec[2] ) / pow(sum_sqrs,3.5)));
    return res;
}


double duz_dy(double* vec, double mu, double J){
    double sum_sqrs = pow(vec[0],2) + pow(vec[1],2) + pow(vec[2],2);
    double res = -mu * ((3*vec[1] * vec[2])/pow(sum_sqrs,2.5)) +
                 J * 1.5 * (((- 8 * pow(vec[0], 2) * vec[1] * vec[2] - 8 * pow(vec[1], 3) * vec[2] + 22 * vec[1] * pow(vec[2], 3)) / pow(sum_sqrs, 4))
                            - (1.0/3.0)* ((15*vec[1]*vec[2] ) / pow(sum_sqrs,3.5)));
    return res;
}


double duz_dz(double* vec, double mu, double J){
    double sum_sqrs = pow(vec[0],2) + pow(vec[1],2) + pow(vec[2],2);
    double res = -mu * ((2* pow(vec[2],2) - pow(vec[0],2) - pow(vec[1], 2))/pow(sum_sqrs,2.5)) +
                 J * 1.5 * (((2*pow(vec[0],4) + 4*pow(vec[0],2)*pow(vec[1],2) - 19*pow(vec[0],2)*pow(vec[2],2)
                 + 2*pow(vec[1],4) - 19*pow(vec[1],2)*pow(vec[2],2) + 9*pow(vec[2],2)) / pow(sum_sqrs,4))
                 - (1.0/3.0)* ((15*pow(vec[2],2) * pow(sum_sqrs,0.5) - 3* pow(sum_sqrs, 1.5)) / pow(sum_sqrs,4)));
    return res;

}

double dux_dmu(double* vec){
    double sum_sqrs = pow(vec[0],2) + pow(vec[1],2) + pow(vec[2],2);
    double res = vec[0] / pow(sum_sqrs, 1.5);
    return res;
}

double duy_dmu(double* vec){
    double sum_sqrs = pow(vec[0],2) + pow(vec[1],2) + pow(vec[2],2);
    double res = vec[1] / pow(sum_sqrs, 1.5);
    return res;
}

double duz_dmu(double* vec){
    double sum_sqrs = pow(vec[0],2) + pow(vec[1],2) + pow(vec[2],2);
    double res = vec[2] / pow(sum_sqrs, 1.5);
    return res;
}

double dux_dJ(double* vec){
    double sum_sqrs = pow(vec[0],2) + pow(vec[1],2) + pow(vec[2],2);
    double res = 1.5 * (( (-5*vec[0]*pow(vec[2],2))/pow(sum_sqrs, 3.5))
                        - (1.0/3.0)* ( -3*vec[0] / pow(sum_sqrs, 2.5) ));
    return res;
}

double duy_dJ(double* vec){
    double sum_sqrs = pow(vec[0],2) + pow(vec[1],2) + pow(vec[2],2);
    double res = 1.5 * (( (-5*vec[1]*pow(vec[2],2))/pow(sum_sqrs, 3.5))
                        - (1.0/3.0)* ( -3*vec[1] / pow(sum_sqrs, 2.5) ));
    return res;
}

double duz_dJ(double* vec){
    double sum_sqrs = pow(vec[0],2) + pow(vec[1],2) + pow(vec[2],2);
    double res = 1.5 * (( (2*vec[2]*pow(sum_sqrs,1.5) - 5* pow(vec[2],3)* pow(sum_sqrs, 0.5))/pow(sum_sqrs, 4))
                        - (1.0/3.0)* ( -3*vec[2] / pow(sum_sqrs, 2.5) ));
    return res;
}

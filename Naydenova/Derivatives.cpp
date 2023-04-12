#include "Derivatives.h"


double dx(double* vec) {
    /*
    double res = -MU_CONST * ((-vec[0] / ((pow(vec[0],2) + pow(vec[1],2) + pow(vec[2],2)) * sqrt(pow(vec[0],2) + pow(vec[1],2) + pow(vec[2],2))))
            + J2 * 1.5 * ( (( -6.0 * pow(vec[2],2)*vec[0]) / pow(pow(vec[0],2) + pow(vec[1],2) + pow(vec[2],2) ,4))
            - (1.0/3.0) * ((-4.0 *vec[0]) / (pow(pow(vec[0],2) + pow(vec[1],2) + pow(vec[2],2),3)))));
    */
    double sum_sqrs = pow(vec[0],2) + pow(vec[1],2) + pow(vec[2],2);
    double res = -MU_CONST * (-vec[0] / pow(sum_sqrs, 1.5)) +
            J2 * 1.5 * (( (-5*vec[0]*pow(vec[2],2))/pow(sum_sqrs, 3.5))
            - (1.0/3.0)* ( -3*vec[0] / pow(sum_sqrs, 2.5) ));
    return res;
}

double dy(double* vec) {
    /*
    double res = -MU_CONST * ((-vec[1] / ((pow(vec[0],2) + pow(vec[1],2) + pow(vec[2],2)) * sqrt(pow(vec[0],2) + pow(vec[1],2) + pow(vec[2],2))))
            + J2 * 1.5 * ( (( -6.0 * pow(vec[2],2) * vec[1] ) / pow(pow(vec[0],2) + pow(vec[1],2) + pow(vec[2],2) ,4))
            - (1.0/3.0) * ((-4.0 *vec[1]) / (pow(pow(vec[0],2) + pow(vec[1],2) + pow(vec[2],2),3)))));
            */
    double sum_sqrs = pow(vec[0],2) + pow(vec[1],2) + pow(vec[2],2);
    double res = -MU_CONST * (-vec[1] / pow(sum_sqrs, 1.5)) +
                 J2 * 1.5 * (( (-5*vec[1]*pow(vec[2],2))/pow(sum_sqrs, 3.5))
                 - (1.0/3.0)* ( -3*vec[1] / pow(sum_sqrs, 2.5) ));

    return res;
}

double dz(double* vec) {
    /*
    double res = -MU_CONST * ((-vec[2] / ((pow(vec[0],2) + pow(vec[1],2) + pow(vec[2],2)) * sqrt(pow(vec[0],2) + pow(vec[1],2) + pow(vec[2],2))))
            + J2 * 1.5 * ( ( (2.0 *pow(vec[0],2)*vec[2] + 2*pow(vec[1],2)*vec[2] - 4 * pow(vec[2], 3)) / pow(pow(vec[0],2) + pow(vec[1],2) + pow(vec[2],2) ,4))
            - (1.0/3.0) * ((-4.0 *vec[2]) / (pow(pow(vec[0],2) + pow(vec[1],2) + pow(vec[2],2),3)))));
    */
    double sum_sqrs = pow(vec[0],2) + pow(vec[1],2) + pow(vec[2],2);
    double res = -MU_CONST * (-vec[2] / pow(sum_sqrs, 1.5)) +
                 J2 * 1.5 * (( (2*vec[2]*pow(sum_sqrs,1.5) - 5* pow(vec[2],3)* pow(sum_sqrs, 0.5))/pow(sum_sqrs, 4))
                 - (1.0/3.0)* ( -3*vec[2] / pow(sum_sqrs, 2.5) ));

    return res;
}

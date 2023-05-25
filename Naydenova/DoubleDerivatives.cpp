#include "DoubleDerivatives.h"



double dux_dx(double* vec, double mu, double J){

    double sum_sqrs = pow(vec[0],2) + pow(vec[1],2) + pow(vec[2],2);
    /*
    double res = -mu * ((2.0* pow(vec[0],2) - pow(vec[1],2) - pow(vec[2], 2))/pow(sum_sqrs,2.5)) +
                 J * 1.5 * (((5.0 *pow(vec[2],2)*(6.0*pow(vec[0],2) - pow(vec[1],2) - pow(vec[2],2))) / pow(sum_sqrs,4.5))
                 - (1.0/3.0)* ((3.0 * (4 * pow(vec[0], 2) - pow(vec[1], 2) -  pow(vec[2], 2))) / pow(sum_sqrs,3.5))); */
    double res = 1.5*J*(35.0*pow(vec[0],2)*pow(vec[2],2)/pow(sum_sqrs,4.5) - 5.0*pow(vec[0],2)/pow(sum_sqrs,3.5) - 5.0*pow(vec[2],2)/pow(sum_sqrs,3.5) + 1.0/pow(sum_sqrs,2.5)) - 3.0*mu*pow(vec[0],2)/pow(sum_sqrs, 2.5) + 1.0*mu/pow(sum_sqrs,1.5);
    return -res;
}

double dux_dy(double* vec, double mu, double J){

    double sum_sqrs = pow(vec[0],2) + pow(vec[1],2) + pow(vec[2],2);
    /*
    double res = -mu * ((3.0*vec[0] * vec[1])/pow(sum_sqrs,2.5)) +
                 J * 1.5 * (((35.0*vec[0]* vec[1] * pow(vec[2],2))/ (pow(sum_sqrs,4.5)))
                 - (1.0/3.0)* ((15.0*vec[0]*vec[1] ) / pow(sum_sqrs,3.5)));*/
    double res = 1.5*J*(35.0*vec[0]*vec[1]*pow(vec[2],2)/pow(sum_sqrs,4.5) - 5.0*vec[0]*vec[1]/pow(sum_sqrs,3.5)) - 3.0*mu*vec[0]*vec[1]/pow(sum_sqrs, 2.5);

    return -res;
}

double dux_dz(double* vec, double mu, double J){

    double sum_sqrs = pow(vec[0],2) + pow(vec[1],2) + pow(vec[2],2);
    /*
    double res = -mu * ((3*vec[0] * vec[2])/pow(sum_sqrs,2.5))  +
                 J * 1.5 * (((-5.0*vec[0]*vec[2]*(2.0 * pow(vec[0],2) + 2.0 * pow(vec[1],2) - 5.0*pow(vec[2],2))) / pow(sum_sqrs,4.5))
                 - (1.0/3.0)* ((15.0*vec[0]*vec[2] ) / pow(sum_sqrs,3.5)));*/
    double res = 1.5*J*(35.0*vec[0]*pow(vec[2],3)/pow(sum_sqrs,4.5) - 15.0*vec[0]*vec[2]/pow(sum_sqrs,3.5)) - 3.0*mu*vec[0]*vec[2]/pow(sum_sqrs, 2.5);

    return -res;
}

double duy_dx(double* vec, double mu, double J){

    double sum_sqrs = pow(vec[0],2) + pow(vec[1],2) + pow(vec[2],2);
    /*
    double res = -mu * ((3.0*vec[0] * vec[1])/pow(sum_sqrs,2.5)) +
                 J * 1.5 * (((35.0 *vec[0]* vec[1] * pow(vec[2],2))/ (pow(sum_sqrs,4.5)))
                 - (1.0/3.0)* ((15.0 *vec[0]*vec[1] ) / pow(sum_sqrs,3.5)));*/
    double res = 1.5*J*(35.0*vec[0]*vec[1]*pow(vec[2],2)/pow(sum_sqrs,4.5) - 5.0*vec[0]*vec[1]/pow(sum_sqrs,3.5)) - 3.0*mu*vec[0]*vec[1]/pow(sum_sqrs, 2.5);

    return -res;
}


double duy_dy(double* vec, double mu, double J){

    double sum_sqrs = pow(vec[0],2) + pow(vec[1],2) + pow(vec[2],2);
    /*
    double res = -mu * ((2.0 * pow(vec[1],2) - pow(vec[0],2) - pow(vec[2], 2))/pow(sum_sqrs,2.5)) +
                 J * 1.5 * (((-5.0 *pow(vec[2],2)* (pow(vec[0],2) -6.0 * pow(vec[1],2) + pow(vec[2],2))) / pow(sum_sqrs,4.5))
                 - (1.0/3.0)* ((-3.0 * (pow(vec[0], 2) - 4.0 * pow(vec[1], 2) + pow(vec[2], 2))) / pow(sum_sqrs,3.5)));*/
    double res = 1.5*J*(35.0*pow(vec[1],2)*pow(vec[2],2)/pow(sum_sqrs,4.5) - 5.0*pow(vec[1],2)/pow(sum_sqrs,3.5) - 5.0*pow(vec[2],2)/pow(sum_sqrs,3.5) + 1.0/pow(sum_sqrs,2.5)) - 3.0*mu*pow(vec[1],2)/pow(sum_sqrs, 2.5) + 1.0*mu/pow(sum_sqrs,1.5);

    return -res;
}


double duy_dz(double* vec, double mu, double J){

    double sum_sqrs = pow(vec[0],2) + pow(vec[1],2) + pow(vec[2],2);
    /*
    double res = -mu * ((3.0 *vec[1] * vec[2])/pow(sum_sqrs,2.5))  +
                 J * 1.5 * (((-5.0 *vec[1]*vec[2]*(2.0 * pow(vec[0],2) + 2.0 * pow(vec[1],2) - 5.0 *pow(vec[2],2))) / pow(sum_sqrs,4.5))
                 - (1.0/3.0)* ((15.0 *vec[1]*vec[2] ) / pow(sum_sqrs,3.5)));*/
   double res = 1.5*J*(35.0*vec[1]*pow(vec[2],3)/pow(sum_sqrs,4.5) - 15.0*vec[1]*vec[2]/pow(sum_sqrs,3.5)) - 3.0*mu*vec[1]*vec[2]/pow(sum_sqrs, 2.5);

    return -res;
}


double duz_dx(double* vec, double mu, double J){

    double sum_sqrs = pow(vec[0],2) + pow(vec[1],2) + pow(vec[2],2);
    /*
    double res = -mu * ((3.0 *vec[0] * vec[2])/pow(sum_sqrs,2.5)) +
                 J * 1.5 * (((-5.0 *vec[0]*vec[2]* (2.0 * pow(vec[0], 2) + 2.0*pow(vec[1],2) - 5.0*pow(vec[2], 2))) / pow(sum_sqrs, 4.5))
                 - (1.0/3.0)* ((15.0*vec[0]*vec[2] ) / pow(sum_sqrs,3.5))); */
    double res = 1.5*J*(35.0*vec[0]*pow(vec[2],3)/pow(sum_sqrs,4.5) - 15.0*vec[0]*vec[2]/pow(sum_sqrs,3.5)) - 3.0*mu*vec[0]*vec[2]/pow(sum_sqrs, 2.5);

    return -res;
}


double duz_dy(double* vec, double mu, double J){

    double sum_sqrs = pow(vec[0],2) + pow(vec[1],2) + pow(vec[2],2);
    /*
       double res = -mu * ((3.0*vec[1] * vec[2])/pow(sum_sqrs,2.5)) +
                     J * 1.5 * (((-5.0*vec[1]*vec[2]* (2.0 * pow(vec[0], 2) + 2.0*pow(vec[1],2) - 5.0*pow(vec[2], 2))) / pow(sum_sqrs, 4.5))
                     - (1.0/3.0)* ((15.0*vec[1]*vec[2] ) / pow(sum_sqrs,3.5)));
                     */
    double res = 1.5*J*(35.0*vec[1]*pow(vec[2],3)/pow(sum_sqrs,4.5) - 15.0*vec[1]*vec[2]/pow(sum_sqrs,3.5)) - 3.0*mu*vec[1]*vec[2]/pow(sum_sqrs, 2.5);
    return -res;
}


double duz_dz(double* vec, double mu, double J){

    double sum_sqrs = pow(vec[0],2) + pow(vec[1],2) + pow(vec[2],2);
/*
    double res = -mu * ((2.0* pow(vec[2],2) - pow(vec[0],2) - pow(vec[1], 2))/pow(sum_sqrs,2.5)) +
                 J * 1.5 * (((2.0 * pow(vec[0], 4) + pow(vec[0], 2)* (4.0*pow(vec[1], 2) - 21.0 * pow(vec[2], 2)) +
                 2*pow(vec[1], 4) - 21.0*pow(vec[1], 2) * pow(vec[2], 2) + 12.0* pow(vec[2], 4))/ pow(sum_sqrs, 4.5)  )
                 - (1.0/3.0)* ((-3.0 * (pow(vec[0], 2) + pow(vec[1], 2) - 4.0 * pow(vec[2], 2))) / pow(sum_sqrs,3.5)));
                 */
    double res = 1.5*J*(35.0*pow(vec[2],4)/pow(sum_sqrs,4.5) - 30.0*pow(vec[2],2)/pow(sum_sqrs,3.5) + 3.0/pow(sum_sqrs,2.5)) - 3.0*mu*pow(vec[2],2)/pow(sum_sqrs, 2.5) + 1.0*mu/pow(sum_sqrs,1.5);
    return -res;

}

double dux_dmu(double* vec){

    double sum_sqrs = pow(vec[0],2) + pow(vec[1],2) + pow(vec[2],2);
    /*
   double res = vec[0] / pow(sum_sqrs, 1.5);
     */
    double res = 1.0*vec[0]/pow(sum_sqrs,1.5);
   return -res;
}

double duy_dmu(double* vec){

    double sum_sqrs = pow(vec[0],2) + pow(vec[1],2) + pow(vec[2],2);
    /*
    double res = vec[1] / pow(sum_sqrs, 1.5);*/
    double res = 1.0*vec[1]/pow(sum_sqrs,1.5);
    return -res;
}

double duz_dmu(double* vec){

    double sum_sqrs = pow(vec[0],2) + pow(vec[1],2) + pow(vec[2],2);
    /*
    double res = vec[2] / pow(sum_sqrs, 1.5);
     */
    double res = 1.0*vec[2]/pow(sum_sqrs,1.5);

    return -res;
}

double dux_dJ(double* vec){

    double sum_sqrs = pow(vec[0],2) + pow(vec[1],2) + pow(vec[2],2);
    /*
   double res = 1.5 * (( (-5.0*vec[0]*pow(vec[2],2))/pow(sum_sqrs, 3.5))
                       - (1.0/3.0)* ( (-3.0 *vec[0]) / pow(sum_sqrs, 2.5) ));
                       */
    double res = -7.5*vec[0]*pow(vec[2],2)/pow(sum_sqrs,3.5) + 1.5*vec[0]/pow(sum_sqrs, 2.5);

    return -res;
}

double duy_dJ(double* vec){

    double sum_sqrs = pow(vec[0],2) + pow(vec[1],2) + pow(vec[2],2);
    /*
   double res = 1.5 * (( (-5.0*vec[1]*pow(vec[2],2))/pow(sum_sqrs, 3.5))
                       - (1.0/3.0)* ( (-3.0*vec[1]) / pow(sum_sqrs, 2.5) ));
                       */
    double res = -7.5*vec[1]*pow(vec[2],2)/pow(sum_sqrs,3.5) + 1.5*vec[1]/pow(sum_sqrs, 2.5);

    return -res;
}

double duz_dJ(double* vec){

   double sum_sqrs = pow(vec[0],2) + pow(vec[1],2) + pow(vec[2],2);
    /*
    double res = 1.5 * (((vec[2]* (2.0 * pow(vec[0],2) + 2.0 * pow(vec[1],2) - 3.0 * pow(vec[2],2)))/pow(sum_sqrs, 3.5))
                        - (1.0/3.0)* ( (-3.0 *vec[2]) / pow(sum_sqrs, 2.5) ));
                        */
    double res = -7.5*pow(vec[2],3)/pow(sum_sqrs,3.5) + 4.5*vec[2]/pow(sum_sqrs, 2.5);
    return -res;
}

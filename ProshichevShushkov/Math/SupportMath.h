#ifndef COMPMATH_AE_FACTORIAL_H
#define COMPMATH_AE_FACTORIAL_H

long factorial(int N){
    if(N < 0)
        return 0;
    if (N == 0)
        return 1;
    else
        return N * fact(N - 1);
}

#endif

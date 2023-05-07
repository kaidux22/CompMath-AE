#ifndef COMPMATH_FACTORIAL
#define COMPMATH_FACTORIAL

double factorial(int N){
    if(N < 0)
        return 0.0;
    if (N == 0)
        return 1.0;
    else
        return (double)N * fact(N - 1);
}

#endif
#include "SupportMath.h"

double Factorial(int N){
	double res = 1;
	for(int i = 2; i <= N; i++){
		res *= (double)i;
	}
	return res;

}

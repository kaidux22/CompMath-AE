#ifndef COMPMATH_CONVERTER
#define COMPMATH_CONVERTER

#include "sofa/sofa.h"

/*
Функция транспонирует матрицу 3х3
*/
void Trans(double (*matrix)[3]) {
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			swap(matrix[i][j], matrix[j][i]);
		}
	}
}

/*
Функция применяет матрицу для смены координат
*/
void BodySpaceFixed(double* vec, double time, double (*rotateMatrix)[3]) {
	return;
}

#endif //COMPMATH_CONVERTER
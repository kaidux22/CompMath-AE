#ifndef COMPMATH_LEASTSQUARES
#define COMPMATH_LEASTSQUARES

#include <iostream>
#include <cmath>

#include "../Math/Matrix.h"

/* 
В прямой задаче использовались первые 4 члена ряда гравитационного потенцила.
Тогда количетсво параметров, на которые нужно наложить шум, следующие:
	- Начальная позиция и скорость первого спустника (6 параметров)
	- Начальная позиция и скорость второго спустника (6 параметров)
	- Коэффициенты C и S
		-- при n = 0 (0 параметров)
		-- при n = 2 (6 параметров)
		-- при n = 3 (8 параметров)
		-- при n = 4 (10 параметров)
	- Коэффициент M 
    - Коэффициенты S при m = 0 равны 0, их не учитываем
Итого: 34 неизвестных параметров
*/
#define UNKNOWN_PARAM 34

/* Метод наименьших квадратов */
class LeastSquare{
public:

    /* Иницилизация начальных параметров */
    LeastSquare(double *measure, int measureCnt);

    /* Один шаг метода Ньютона-Гаусса */
    void Iteration(int staps = 1);

    ~LeastSquare();
private:
    int mMeasureCount;

    double *mMeasure;
    double *mStates;
    Matrix<double> *mParams;

    // Невязки
    Matrix<double> *mResiduals;
    Matrix<double> *mMatrixA;
};

#endif //COMPMATH_LEASTSQUARES
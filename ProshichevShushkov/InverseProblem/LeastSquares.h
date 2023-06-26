#ifndef COMPMATH_LEASTSQUARES
#define COMPMATH_LEASTSQUARES

#include <iostream>
#include <cmath>
#include <cstdlib>
#include <fstream>

#include "../Math/Matrix.h"
#include "../NumIntegrate/IntegrateFunc.h"
#include "DistanceOrbits.h"

/* Юлианская дата */
#define JD 2451545.0

/* Шаг в 1 минуту */
#define STEP 30.0 //с

/* G * m */
#define GM 398600.4415 // км^3 / с^2

/* радиус Земли + 500 км */
#define START_POINT 6878.0 //км

 /* угол наклона орбиты */
#define ANGLE 11.0

/* Метод наименьших квадратов */
class LeastSquare{
public:

    /* Иницилизация начальных параметров */
    LeastSquare(double *measure, int measureCnt);

    /* Один шаг метода Ньютона-Гаусса */
    void Iteration(int staps = 1);
private:
    /* подсчёт матрицы dg/dx */
    Matrix<double> *MatrixdGdX();

    /* Решение СЛАУ методом Холецкого */
    Matrix<double> *CholeskyDecomposition(Matrix<double> *MatrixA, Matrix<double> *Vectorb);

    int mMeasureCount; // количество измерений

<<<<<<< HEAD
    double *mMeasure;
    double *mVec;
    Matrix<double> *mStates;
    Matrix<double> *mParams;

    // Невязки
    Matrix<double> *mResiduals;
    Matrix<double> *mMatrixA;
    Matrix<double> *mTruth;
    char* mSymb[UNKNOWN_PARAM] = {"C02", "C12", "C22", 
                                  "C03", "C13", "C23", "C33", 
                                  "C04", "C14", "C24", "C34", "C44", 
                                  "S12", "S22",
                                  "S13", "S23", "S33", 
                                  "S14", "S24", "S34", "S44"};
=======
    double *mMeasure; // измерений
    double *mVec; // расширенный вектор
    Matrix<double> *mStates; // dx/dp
    Matrix<double> *mParams; // восстанавливаемые параметры
    Matrix<double> *mResiduals; //невязки
    Matrix<double> *mMatrixA; //матрица A
    Matrix<double> *mTruth; //истинные значения
    char* mSymb[UNKNOWN_PARAM] = {"C13", "C23", "C33", "C14", "C24", "C34", "C44", "S13", "S23", "S33", "S14", "S24", "S34", "S44"};
>>>>>>> d9a997fa8a9e59c18ef81a71abaa4a30aea08804
};

#endif //COMPMATH_LEASTSQUARES
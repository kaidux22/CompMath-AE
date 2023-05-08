#include "LeastSquares.h"

LeastSquare::LeastSquare(double *measure, int measureCnt){
    mMeasure = measure;
    mMeasureCount = measureCnt;

    //первый индекс - m, второй индекс - n
	//С[0][n] = Jn
	double Cmn[5][5] = { { -1.0, 0.0, 0.1082635854e-2, -0.2532435346e-5, -0.1619331205e-5 },
					   {0.0, 0.0, -0.3504890360e-9, 0.2192798802e-5, -0.5087253036e-6},
					   {0.0, 0.0, 0.1574536043e-5, 0.3090160446e-6, 0.7841223074e-7},
					   {0.0, 0.0, 0.0, 0.1005588574e-6, 0.5921574319e-7},
					   {0.0, 0.0, 0.0, 0.0, -0.3982395740e-8} };
	// первый индекс - m, второй индекс - n
	double Smn[5][5] = { {0.0, 0.0, 0.0, 0.0, 0.0},
					   {0.0, 0.0, 0.1635406077e-8, 0.2680118938e-6, -0.4494599352e-6},
					   {0.0, 0.0, -0.9038680729e-6, -0.2114023978e-6, 0.1481554569e-6},
					   {0.0, 0.0, 0.0, 0.1972013239e-6, -0.1201129183e-7},
					   {0.0, 0.0, 0.0, 0.0, 0.6525605810e-8} };

    double mNoise[34] = {0, 0, 0, 0, 0, 0, 0,
						0, 0, 0, 0, 0, 0, 0,
						0, 0, 0, 0, 0, 0, 0,
						0, 0, 0, 0, 0, 0, 0,
						0, 0, 0, 0, 0, 0};

    // Вектор состояния для одного спутника (начальных 6 параметров | единичная матрица 6х6 | нулевые столбцы для оставшихся 25ти коэффициентов)
    mStates = new double[6 + 6 * 6 + 22 * 6];
    mParams = new Matrix<double>(UNKNOWN_PARAM, 1);
    mResiduals = new Matrix<double>(mMeasureCount, 1);
    mMatrixA = new Matrix<double>(mMeasureCount, UNKNOWN_PARAM);

    // параметры первого спутника
    mParams->Set(0, 0, 6878.0 + mNoise[0]), mParams->Set(1, 0, mNoise[1]), mParams->Set(2, 0, mNoise[2]);
    mParams->Set(3, 0, mNoise[3]), mParams->Set(4, 0, sqrt(398600.4415 / 6878.0) + mNoise[4]), mParams->Set(5, 0, mNoise[5]);

    //параметры второго спутника
    mParams->Set(6, 0, 1472.62 + mNoise[6]), mParams->Set(7, 0, -6718.5 + mNoise[7]), mParams->Set(8, 0, -0.148523 + mNoise[8]);
    mParams->Set(9, 0, 7.43608 + mNoise[9]), mParams->Set(10, 0, 1.63027 + mNoise[10]), mParams->Set(11, 0, 0.000242242 + mNoise[11]);
    
    //Нахождение параметра массы
    mParams->Set(12, 0, 398600.4415 + mNoise[12]);

    int cnt = 13;
    //Jn
    for(int n = 2; n < 5; n++){
        mParams->Set(cnt, 0, Cmn[0][n] + mNoise[cnt]);
        cnt++;
    }

    for(int n = 2; n < 5; n++){
        for(int m = 1; m <= n; m++){
            mParams->Set(cnt, 0, Cmn[m][n] + mNoise[cnt]);
            mParams->Set(cnt + 1, 0, Smn[m][n] + mNoise[cnt + 1]);
            cnt += 2;
        }
    }
}

void LeastSquare::Iteration(int steps){
    for(int step = 0 ; step < steps; step++){
        for(int i = 6; i < 174; i++)
            mStates[i] = 0;

        for(int i = 0; i < 6; i++){
            mStates[i] = mParams->Get(i, 0);
            mStates[6 + i * 7] = 1;
        }

        ConditionVectorIntegrate(JD, STEP, 174, mStates, mParams);
        
    }
}

LeastSquare::~LeastSquare(){
    delete[] mStates;
    delete mParams;
    delete mResiduals;
    delete mMatrixA;
}


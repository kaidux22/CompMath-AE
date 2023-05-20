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

    double mNoise[34];
    srand(time(0));

    for(int i = 0; i < 34; i++){
        mNoise[i] = (double)(rand() % (int)2e9 - 1e9) / 1e9 / 2e3 + 1;
    }

    // Вектор состояния для одного спутника (начальных 12 параметров | единичная матрица 12х12 | нулевые столбцы для оставшихся 22ух коэффициентов)
    mVec = new double[12 + 12 * UNKNOWN_PARAM];
    mStates = new Matrix<double>(12, UNKNOWN_PARAM);
    mParams = new Matrix<double>(UNKNOWN_PARAM, 1);
    mResiduals = new Matrix<double>(mMeasureCount, 1);
    mMatrixA = new Matrix<double>(mMeasureCount, UNKNOWN_PARAM);
    mTruth = new Matrix<double>(UNKNOWN_PARAM, 1);

    // координаты первого и второго спутников
    mParams->Set(0, 0, 1248.77), mParams->Set(1, 0, -6763.69), mParams->Set(2, 0, -0.155766);
    mParams->Set(3, 0, 1472.62), mParams->Set(4, 0, -6718.5), mParams->Set(5, 0, -0.148523);

    //скорости первого и второго спутников
    mParams->Set(6, 0, 7.48616), mParams->Set(7, 0, 1.38216), mParams->Set(8, 0, 0.00024043);
    mParams->Set(9, 0, 7.43608), mParams->Set(10, 0, 1.63027), mParams->Set(11, 0, 0.000242242);
    
    //Нахождение параметра массы
    //mParams->Set(12, 0, 398600.4415 + mNoise[12]);


    int cnt = 12;
    //Cmn
    for(int n = 2; n < 4; n++){
        for(int m = 0; m <= n; m++){
            mParams->Set(cnt, 0, Cmn[m][n] * mNoise[cnt]);
            cnt++;
        }
    }

    //Smn
    for(int n = 2; n <= 4; n++){
        for(int m = 1; m <= n; m++){
            mParams->Set(cnt, 0, Smn[m][n] * mNoise[cnt]);
            cnt++;
        }
    }
    
    
    for(int i = 0; i < UNKNOWN_PARAM; i++){
        mTruth->Set(i, 0, mParams->Get(i, 0));
        mParams->Set(i, 0, mParams->Get(i, 0) * mNoise[i]);
    }

    for(int i = 0; i < UNKNOWN_PARAM; i++)
        cout << mParams->Get(i, 0) << "\t\t" << mTruth->Get(i, 0) << endl;
    cout << endl;

}

void LeastSquare::Iteration(int steps){


    for(int step = 0 ; step < steps; step++){
        for(int i = 0; i < 12; i++){
            for(int j = 0; j < UNKNOWN_PARAM; j++){
                if(i == j){
                    mStates->Set(i, j, 1);
                    continue;
                }
                mStates->Set(i, j, 0);
            }
        }

        for(int i = 0; i < 12; i++){
            mVec[i] = mParams->Get(i, 0);
        }
        
        //проверить переворот матрицы в список по столбцам
        for(int i = 12; i < 12 + 12 * UNKNOWN_PARAM; i++){
            mVec[i] = mStates->TransToVector()[i - 12];
        }

        double **orbits = ConditionVectorIntegrate(JD, STEP, 12 + 12 * UNKNOWN_PARAM, mVec, mParams);

        /*
        for(int t = 0; t < mMeasureCount; t++){
            for(int i = 13; i < 12 + 12 * 34 + 1; i++){
               cout << orbits[t][i] << " ";
            }
            cout << endl;
        }
        cout << endl;
        assert(false);
        */

        double *distance = OrbitDistance(orbits, mMeasureCount);

        /*
        for(int i = 0; i < mMeasureCount; i++){
            cout << distance[2 * i] - JD << " " << distance[2 * i + 1] << endl;
        }
        assert(false);
        */

        for(int i = 0; i < mMeasureCount; i++){
            for(int j = 0; j < 12; j++){
                mVec[j] = orbits[i][j + 1];
            }

            for(int j = 0; j < UNKNOWN_PARAM * 12; j++){
                mStates->TransToVector()[j] = orbits[i][13 + j];
            }

            Matrix<double> *dGdX = MatrixdGdX();
            Matrix<double> prod = (*dGdX * *mStates);                                         
            double *res = prod.TransToVector();

            for(int j = 0; j < UNKNOWN_PARAM; j++){
                mMatrixA->Set(i, j, -res[j]);
            }

            //разобраться со слау
            mResiduals->Set(i, 0, mMeasure[2 * i + 1] - distance[2 * i + 1]);        
        }

        Matrix<double> MatrixAtA(UNKNOWN_PARAM, UNKNOWN_PARAM);
        MatrixAtA = (mMatrixA->Transposition() * *mMatrixA);
 

        Matrix<double> MatrixArb(UNKNOWN_PARAM, 1);
        MatrixArb = mMatrixA->Transposition() * *mResiduals;

        Matrix<double> *Vectorx = CholeskyDecomposition(&MatrixAtA, &MatrixArb);

        for(int i = 0; i < UNKNOWN_PARAM; i++){
            mParams->Set(i, 0, mParams->Get(i, 0) - Vectorx->Get(i, 0));
        }

        for(int i = 0; i < UNKNOWN_PARAM; i++)
            cout << mParams->Get(i, 0) << "\t\t" << mTruth->Get(i, 0) << endl;
        cout << endl;
        //mTruth->Print();

    }
}

Matrix<double>* LeastSquare::MatrixdGdX(){
    Matrix<double> *dGdX = new Matrix<double>(1, 12);

    for(int i = 0; i < 3; i++){
        dGdX->Set(0, i, (mVec[i] - mVec[i + 3]) / sqrt(pow(mVec[0] - mVec[3], 2) + pow(mVec[1] - mVec[4], 2) + pow(mVec[2] - mVec[5], 2)));
        dGdX->Set(0, i + 3, (mVec[i + 3] - mVec[i]) / sqrt(pow(mVec[0] - mVec[3], 2) + pow(mVec[1] - mVec[4], 2) + pow(mVec[2] - mVec[5], 2)));
    }
    return dGdX;
}

Matrix<double>* LeastSquare::CholeskyDecomposition(Matrix<double> *MatrixA, Matrix<double> *Vectorb){
    assert(MatrixA->RowsCount() == MatrixA->ColumnCount());
    assert(MatrixA->RowsCount() == Vectorb->RowsCount());
    Matrix<double> *MatrixL = new Matrix<double>(MatrixA->RowsCount(), MatrixA->ColumnCount());

     for (int i=0; i < MatrixA->RowsCount(); i++){
        for (int j = 0; j < (i + 1); j++){
            double res = 0;
            for (int k = 0; k < j; k++) {
                res += MatrixL->Get(i, k) * MatrixL->Get(j, k);
            }
            if (i == j) {
                MatrixL->Set(i, j, sqrt(MatrixA->Get(i, i) - res));
            } else {
                if (MatrixL->Get(j, j) == 0){
                    MatrixL->Set(i, j, 0);
                }
                else {
                    MatrixL->Set(i, j, (1.0 / MatrixL->Get(j, j) * (MatrixA->Get(i, j) - res)));
                }
            }
        }
    }

    Matrix<double> *Vectorx = new Matrix<double>(Vectorb->RowsCount(), 1);
    Matrix<double> *Vectory = new Matrix<double>(Vectorb->RowsCount(), 1);

    //  L*y=b
    for (int i = 0; i < Vectorb->RowsCount(); i++){
        double res = 0;
        for (int j = 0; j < i; j++){
            res += MatrixL->Get(i, j) * Vectory->Get(j, 0);
        }
        if (MatrixL->Get(i, i) == 0){
            Vectory->Set(i, 0, 0);
        }
        else {
            Vectory->Set(i, 0, (1.0 / MatrixL->Get(i, i)) * (Vectorb->Get(i, 0) - res));
        }
    }

    //  L^t*x=y
    for (int i = Vectorb->RowsCount() - 1; i >= 0; i--){
        double res = 0;
        for (int j = i+1; j < Vectorb->RowsCount(); j++){
            res += MatrixL->Get(j, i) * Vectorx->Get(j, 0);
        }

        Vectorx->Set(i, 0, (1.0 / MatrixL->Get(i, i)) * (Vectory->Get(i, 0) - res));

    }

    for (int i = Vectorb->RowsCount() - 1; i >= 0; i--){
        double res = 0;
        for (int j = i+1; j < Vectorb->RowsCount(); j++){
            res += MatrixL->Get(j, i) * Vectorx->Get(j, 0);
        }
        if (MatrixL->Get(i, i) == 0){
            Vectorx->Set(i, 0, 0);
        }
        else {
            Vectorx->Set(i, 0, (1.0 / MatrixL->Get(i, i) * (Vectory->Get(i, 0) - res)));
        }
 
 
    }

    delete Vectory;
    delete MatrixL;
    return Vectorx;
}

LeastSquare::~LeastSquare(){
    /*
    delete[] mVec;
    delete mStates;
    delete mParams;
    delete mResiduals;
    delete mMatrixA;
    */
}




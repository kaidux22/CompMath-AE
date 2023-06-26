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

    double mNoise[UNKNOWN_PARAM];
    srand(time(0));

<<<<<<< HEAD
    for(int i = 0; i < UNKNOWN_PARAM; i++){
        mNoise[i] = 1.0; //(double)(rand() % (int)2e9 - 1e9) / 1e9 / 2e2 + 1;
=======
    //генерируем шум
    for(int i = 0; i < UNKNOWN_PARAM; i++){
        mNoise[i] = (rand() % (int)2e5 - 1e5) / 1e8 + 1;
>>>>>>> d9a997fa8a9e59c18ef81a71abaa4a30aea08804
    }

    // начальных 12 параметров | нулевые столбцы для восстанавливаемых 14ти коэффициентов
    mVec = new double[12 + 12 * UNKNOWN_PARAM];
    mStates = new Matrix<double>(12, UNKNOWN_PARAM);
    mParams = new Matrix<double>(UNKNOWN_PARAM, 1);
    mResiduals = new Matrix<double>(mMeasureCount, 1);
    mMatrixA = new Matrix<double>(mMeasureCount, UNKNOWN_PARAM);
    mTruth = new Matrix<double>(UNKNOWN_PARAM, 1);
<<<<<<< HEAD

    int cnt = 0;

    //Cmn
    for(int n = 2; n <= 4; n++){
        for(int m = 0; m <= n; m++){
=======
    

    int cnt = 0;

    for(int n = 3; n <= 4; n++){  
        for(int m = 0; m <= n; m++){
            if(m == 0)
                continue;
>>>>>>> d9a997fa8a9e59c18ef81a71abaa4a30aea08804
            mParams->Set(cnt, 0, Cmn[m][n]);
            cnt++;
        }
    }


    for(int n = 3; n <= 4; n++){
        for(int m = 1; m <= n; m++){
            mParams->Set(cnt, 0, Smn[m][n]);
            cnt++;
        }
    }
<<<<<<< HEAD
=======

>>>>>>> d9a997fa8a9e59c18ef81a71abaa4a30aea08804
    
    for(int i = 0; i < UNKNOWN_PARAM; i++){
        mTruth->Set(i, 0, mParams->Get(i, 0));
        mParams->Set(i, 0, mParams->Get(i, 0) * mNoise[i]);
    }
<<<<<<< HEAD
    cout << "\tCurrent Value\t" << "True Value\t" << "Difference\n";
    for(int i = 0; i < UNKNOWN_PARAM; i++){
        cout << mSymb[i] << "\t" << mParams->Get(i, 0) << "\t" << mTruth->Get(i, 0) << "\t" << mParams->Get(i, 0) - mTruth->Get(i, 0) << endl;
    }
    cout << endl;
=======

    cout << "\t" << "Offset value" << "\t" << "True value" << "\t" << "Difference" << endl;
>>>>>>> d9a997fa8a9e59c18ef81a71abaa4a30aea08804
}

void LeastSquare::Iteration(int steps){
    double **orbits;

    for(int step = 0 ; step < steps; step++){
        for(int i = 0; i < 12; i++){
            for(int j = 0; j < UNKNOWN_PARAM; j++){
                mStates->Set(i, j, 0);
            }
        }

<<<<<<< HEAD
       //начальное положение в НСК первого спутника
        mVec[0] = 1248.77, mVec[1] = -117.887, mVec[2] = -6762.66;

        //начальное положение в НСК второго спутника
        mVec[3] = 1472.62, mVec[4] = -117.105, mVec[5] = -6717.48;

        // начальная скорость в НСК первого спутника
        mVec[6] = 7.48616, mVec[7] = 0.0238816, mVec[8] = 1.38195;

        //начальная скорость в НСК второго спутника
        mVec[9] = 7.43608, mVec[10] = 0.0282099, mVec[11] = 1.63003;


        //проверить переворот матрицы в список по столбцам
=======
         // координаты первого и второго спутников
        mVec[0] = 1248.77, mVec[1] = -6763.69, mVec[2] = -0.155766;
        mVec[3] = 1472.62, mVec[4] = -6718.5, mVec[5] = -0.148523;

        //скорости первого и второго спутников
        mVec[6] = 7.48616, mVec[7] = 1.38216, mVec[8] = 0.00024043;
        mVec[9] = 7.43608, mVec[10] = 1.63027, mVec[11] = 0.000242242;

        double angle = ANGLE * M_PI / 180.0;

        double rotateMatrix[3][3] = {{1.0, 0.0, 0.0}, {0.0, cos(angle), -sin(angle)}, {0, sin(angle), cos(angle)}};

	    for(int i = 0; i < 4; i++){
	    	changeCoords(rotateMatrix, mVec, 3 * i);
	    }

>>>>>>> d9a997fa8a9e59c18ef81a71abaa4a30aea08804
        for(int i = 12; i < 12 + 12 * UNKNOWN_PARAM; i++){
            mVec[i] = mStates->TransToVector()[i - 12];
        }

        orbits = ConditionVectorIntegrate(JD, STEP, 12 + 12 * UNKNOWN_PARAM, mVec, mParams);

        double *distance = OrbitDistance(orbits, mMeasureCount);

        //составление матрицы A и невязок
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

            mResiduals->Set(i, 0, mMeasure[2 * i + 1] - distance[2 * i + 1]);        
        }

        //подготовка СЛАУ
        Matrix<double> MatrixAtA(UNKNOWN_PARAM, UNKNOWN_PARAM);
        MatrixAtA = (mMatrixA->Transposition() * *mMatrixA);
        Matrix<double> MatrixArb(UNKNOWN_PARAM, 1);
        MatrixArb = mMatrixA->Transposition() * *mResiduals;

        //решение СЛАУ
        Matrix<double> *Vectorx = CholeskyDecomposition(&MatrixAtA, &MatrixArb);

        //переход к вновь восстановленным параметрам
        for(int i = 0; i < UNKNOWN_PARAM; i++){
            mParams->Set(i, 0, mParams->Get(i, 0) - Vectorx->Get(i, 0));
        }

<<<<<<< HEAD
        for(int i = 0; i < UNKNOWN_PARAM; i++){
        cout << mSymb[i] << "\t" << mParams->Get(i, 0) << "\t" << mTruth->Get(i, 0) << "\t" << mParams->Get(i, 0) - mTruth->Get(i, 0) << endl;
        }
=======
        for(int i = 0; i < UNKNOWN_PARAM; i++)
            cout << mSymb[i] << "\t" <<  mParams->Get(i, 0) << "\t" << mTruth->Get(i, 0) << "\t" << mParams->Get(i, 0) - mTruth->Get(i, 0) << endl;
>>>>>>> d9a997fa8a9e59c18ef81a71abaa4a30aea08804
        cout << endl;

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

    //составление матрицы L
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




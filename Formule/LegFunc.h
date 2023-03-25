#ifndef COMPMATH_LEGFUN
#define COMPMATH_LEGFUN

#include <cmath>
#include <cassert>
#include <iostream>

using namespace std;

class LegFunc {
public:
	LegFunc(int n, int m, double arg) {
		mOrd = n;
		mPwr = m;
		mFuncs = new double * [mPwr + 1];
		mFuncs[0] = LegPol(arg, mOrd);

		for (int idx = 1; idx <= mPwr; idx++) {
			mFuncs[idx] = new double[mOrd + 1];
			mFuncs[idx][0] = 0;
		}

		for (int pwr = 1; pwr <= mPwr; pwr++) {
			for (int ord = 1; ord <= mOrd; ord++) {
				if (pwr == ord) {
					mFuncs[pwr][ord] = (double)(2 * ord - 1) * sqrt(1 - arg * arg) * mFuncs[ord - 1][ord - 1];
				}
				else {
					if (ord - 1 < pwr) {
						mFuncs[pwr][ord] = 0;
					}
					else if (ord - 2 < pwr) {
						mFuncs[pwr][ord] = (double)(2 * ord - 1) * arg * mFuncs[pwr][ord - 1] / (double)(ord - pwr);
					}
					else {
						mFuncs[pwr][ord] = ((double)(2 * ord - 1) * arg * mFuncs[pwr][ord - 1] - (double)(ord - 1 + pwr) * mFuncs[pwr][ord - 2]) / ((double)(ord - pwr));
					}
				}
			}
		}

	}

	void PrintMaxtrix() {
		for (int pwr = 0; pwr <= mPwr; pwr++) {
			for (int ord = 0; ord <= mOrd; ord++) {
				cout << mFuncs[pwr][ord] << " ";
			}
			cout << endl;
		}
	}

	double ExtractValue(int n, int m) {
		assert(n <= mOrd && m <= mPwr);
		return mFuncs[m][n];
	}

	~LegFunc() {
		for (int idx = 0; idx <= mPwr; idx++) {
			delete[] mFuncs[idx];
		}
		delete[] mFuncs;
	}

private:

	double* LegPol(double arg, int ord) {
		double* pols = new double[ord + 1];
		pols[0] = 1;

		if (ord == 0)
			return pols;

		pols[1] = arg;

		for (int i = 1; i < ord; i++) {
			pols[i + 1] = ((double)(2 * i + 1) * arg * pols[i] - (double)(i)*pols[i - 1]) / (double)(i + 1);
		}
		return pols;
	}

	double** mFuncs;
	int mOrd, mPwr;
};

#endif //COMPMATH_LEGFUN

#include <cmath>

class LegFunc {
public:
	LegFunc(int n, int m, double arg) {
		mOrd = n;
		mPwr = m;
		mFuncs = new double* [mPwr + 1];
		mFuncs[0] = LegPol(arg, mOrd);
		
		for (int idx = 1; idx <= mPwr; idx++) {
			mFuncs[idx] = new double[mOrd + 1];
			mFuncs[idx][0] = 0;
		}

		if (!mOrd || !mPwr) {
			return;
		}
		mFuncs[1][1] = sqrt(1 - arg * arg);

		for (int ord = 2; ord <= mOrd; ord++) {
			mFuncs[1][ord] = ((double)(2 * mOrd + 3) * arg * mFuncs[1][ord - 1] - (double)(n + m + 1) * mFuncs[1][ord - 2]) / ((double)(n + 2 - m));
		}

		for (int pwr = 2; pwr <= mPwr; pwr++) {
			for (int ord = 1; ord <= mOrd; ord++) {
				mFuncs[pwr][ord] = ((double)(2 * m + 2) * arg * mFuncs[pwr - 1][ord]) / (sqrt(1 - arg * arg)) - (double)((n - m) * (n + m + 1)) * mFuncs[pwr - 2][ord];
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

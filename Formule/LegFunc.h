#include <cmath>
#include <cassert>

typedef long double ld;

using namespace std;

//класс хранит полиномы Лежандра всех степеней и порядков до m и n
class LegFunc {
public:
	//создаём таблицу значений полиномов 
	LegFunc(int n, int m, ld arg) {
		mOrd = n; //количество столбцов/порядок полинома
		mPwr = m; //количество строк/степень полинома
		mFuncs = new ld* [mPwr + 1];
		mFuncs[0] = LegPol(arg, mOrd); //m = 0, поэтому первая строка - обычные полиномы Лежандра
		
		for (int idx = 1; idx <= mPwr; idx++) {
			mFuncs[idx] = new ld[mOrd + 1];
			mFuncs[idx][0] = 0;
		}
		
		//заполняем таблицу по реккурентным формулам
		for (int pwr = 1; pwr <= mPwr; pwr++) {
			for (int ord = 1; ord <= mOrd; ord++) {
				if (pwr == ord) {
					mFuncs[pwr][ord] = (ld)(2 * ord - 1) * sqrt(1 - arg * arg) * mFuncs[ord - 1][ord - 1];
				}
				else {
					if (ord - 1 < pwr) {
						mFuncs[pwr][ord] = 0;
					}
					else if (ord - 2 < pwr) {
						mFuncs[pwr][ord] = (ld)(2 * ord - 1) * arg * mFuncs[pwr][ord - 1] / (ld)(ord - pwr);
					}
					else {
						mFuncs[pwr][ord] = ((ld)(2 * ord - 1) * arg * mFuncs[pwr][ord - 1] - (ld)(ord - 1 + pwr) * mFuncs[pwr][ord - 2]) / ((ld)(ord - pwr));
					}
				}
			}
		}

	}
	
	//выводим всю таблицу значений
	void PrintMaxtrix() {
		for (int pwr = 0; pwr <= mPwr; pwr++) {
			for (int ord = 0; ord <= mOrd; ord++) {
				cout << mFuncs[pwr][ord] << " ";
			}
			cout << endl;
		}
	}

	//выводим значение полинома при  конкретных n и m
	ld ExtractValue(int n, int m) {
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
	//функция создаёт список из полиномов Лежандра до порядка n
	ld* LegPol(ld arg, int ord) {
		ld* pols = new ld[ord + 1];
		pols[0] = 1;

		if (ord == 0)
			return pols;

		pols[1] = arg;

		for (int i = 1; i < ord; i++) {
			pols[i + 1] = ((ld)(2 * i + 1) * arg * pols[i] - (ld)(i)*pols[i - 1]) / (ld)(i + 1);
		}
		return pols;
	}

	ld** mFuncs;
	int mOrd, mPwr;
};

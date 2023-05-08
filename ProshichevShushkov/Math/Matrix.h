#ifndef COMPMATH_MATRIX
#define COMPMATH_MATRIX

#include <cassert>
#include <iostream>

using namespace std;

/* тип данных: матрицы, хранящая Type-элементы */
template <typename Type>
class Matrix {
public:

	/* конструктор, задающий размеры матрицы */
	Matrix(int n, int m) {
		mRows = n;
		mColumn = m;
		mMatrix = new Type[n * m];
		
		for (int i = 0; i < n * m; i++)
			mMatrix[i] = (Type)0;
	}

	/* получить элемент по строке и столбцу */
	Type Get(int i, int j) const{
		assert(i < mRows && j < mColumn);
		return mMatrix[j * mColumn + i];
	}

	/* записать элемент по строке и столбцу */
	void Set(int i, int j, Type value) {
		assert(i < mRows && j < mColumn);
		mMatrix[j * mColumn + i] = value;
	}

	/* транспонирование матрицы */
	 void Transposition() {
		for (int i = 0; i < mRows; i++) {
			for (int j = i; j < mColumn; j++) {
				swap(mMatrix[j * mColumn + i], mMatrix[i * mColumn + j]);
			}
		}
	}

	/* вывод матрицы */
	void Print() {
		for (int i = 0; i < mRows; i++) {
			for (int j = 0; j < mColumn; j++) {
				cout << mMatrix[j * mColumn + i] << " ";
			}
			cout << endl;
		}
		cout << endl;
	}

	/* произведение матриц */
	Matrix<Type> operator *(const Matrix<Type>& other) {
		assert(mColumn == other.mRows);

		Matrix<Type> newMatrix(mRows, other.mColumn);

		for (int i = 0; i < mRows; i++) {
			for (int j = 0; j < other.mColumn; j++) {

				newMatrix.Set(i, j, 0);
				for (int k = 0; k < mColumn; k++) {
					newMatrix.Set(i, j, newMatrix.Get(i, j) + Get(i, k) * other.Get(k, j));
				}
			}
		}

		return newMatrix;
	}

	/* сумма матриц */
	Matrix<Type> operator +(const Matrix<Type> &other) {
		assert(mRows == other.mRows && mColumn == other.mColumn);

		Matrix<Type> newMatrix(mRows, mColumn);
		
		for (int i = 0; i < mRows; i++) {
			for (int j = 0; j < mColumn; j++) {
				newMatrix.Set(i, j, Get(i, j) + other.Get(i, j));
			}
		}
		return newMatrix;

	}

	/* домножение на скаляр */
	void Product(double num) {

		for (int i = 0; i < mRows; i++) {
			for (int j = 0; j < mColumn; j++) {
				mMatrix[j * mColumn + i] *= num;
			}
		}
	}
		
	/* перегрузка оператора присваивания */
	Matrix<Type>& operator=(Matrix<Type>& matrix) {
		assert(mRows == matrix.mRows && mColumn == matrix.mColumn);

		for (int i = 0; i < mRows; i++) {
			for (int j = 0; j < mColumn; j++) {
				mMatrix[j * mColumn + i] = matrix.Get(i, j);
			}
		}
		return *this;
	}

	/* извлечение строки */
	Matrix<Type> ExtractRow(int rowNum) {
		Matrix<Type> vector(1, mColumn);

		for (int i = 0; i < mColumn; i++) {
			vector.Set(i, 0, Get(rowNum, i));
		}

		return vector;
	}

	/* извлечение столбца */
	Matrix<Type> ExctractColunm(int columnNum) {
		Matrix<Type> vector(1, mRows);

		for (int i = 0; i < mRows; i++) {
			vector.Set(i, 0, Get(columnNum, i));
		}
		
		return vector;
	}

	~Matrix() {
		delete[] mMatrix;
	}

private:
	Type* mMatrix;
	int mRows, mColumn;
};

#endif //COMPMATH_MATRIX

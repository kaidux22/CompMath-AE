#ifndef COMPMATH_COMPLEXNUMS
#define COMPMATH_COMPLEXNUMS

#include <iostream>

using namespace std;

//комплексные числа
class ComplexNum {
public:
	//задаём a+bi, как ComplexNum(a, b)
	ComplexNum(double real, double img): mReal(real), mImg(img) {}

	//вещественная часть
	double Real() const { return mReal; }

	//мнимая часть
	double Img() const { return mImg; }

	//сумма комплексных чисел
	ComplexNum operator +(const ComplexNum& other) const {
		return ComplexNum(mReal + other.Real(), mImg + other.Img());
	}

	//разность комплексных чисел
	ComplexNum operator -(const ComplexNum& other) const {
		return ComplexNum(mReal - other.Real(), mImg - other.Img());
	}

	//произведение комплексных числе
	ComplexNum operator *(const ComplexNum& other) {
		return ComplexNum(mReal * other.Real() - mImg * other.Img(), mReal * other.Img() + mImg * other.Real());
	}

	//вывод комплексных чисел в стандартной форме
	friend std::ostream& operator <<(std::ostream& out, const ComplexNum& num) {
		out << num.Real() << "+" << num.Img() << "i";
		return out;
	}

private:
	double mReal, mImg;
};

#endif //COMPMATH_COMPLEXNUMS

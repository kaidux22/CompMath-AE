#ifndef COMPMATH_COMPLEXNUMS
#define COMPMATH_COMPLEXNUMS

#include <iostream>
#include <cmath>

using namespace std;

//комплексные числа
class ComplexNum {
public:
	//задаём a+bi, как ComplexNum(a, b)
	ComplexNum(double real, double img) : mReal(real), mImg(img) {}

	//преобразования double в ComplexNum
	explicit ComplexNum(double realNum) : ComplexNum(realNum, 0) {}

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
	
	//модуль комплексного числа
	double Module() {
		return sqrt(mReal * mReal + mImg * mImg);
	}

	//агрумент комплексного числа
	double Arg() {
		return atan(mImg / mReal);
	}
	
	//Формула Муавра
	ComplexNum Pow(double k) {
		return ((ComplexNum)pow(Module(), k)) * ComplexNum(cos(Arg() * k), sin(Arg() * k));
	}

	friend std::ostream& operator <<(std::ostream& out, const ComplexNum& num) {
		out << num.Real() << "+" << num.Img() << "i";
		return out;
	}



private:
	double mReal, mImg;
};

#endif //COMPMATH_COMPLEXNUMS

#include "main.h"

using namespace std;

int main()
{	double JD_start = JD;
	double *vec = new double[6];
	double rotateMatrix[3][3];
	//начальное положение в ЗСК
	vec[0] = START_POINT, vec[1] = 0, vec[2] = 0;

	// начальная скорость в ЗСК
	vec[3] = 0, vec[4] = sqrt(GM / START_POINT), vec[5] = 0;

	//перевод в НСК начальных параметров
	iauC2t06a(JD_start + (37.0 + 32.184) / 86400.0, 0, JD_start, 0, 0, 0, rotateMatrix);
	Transposition(rotateMatrix);
	changeCoords(rotateMatrix, vec, 0); //перевод начальных координат в НСК
	changeCoords(rotateMatrix, vec, 3); //перевод проекций скоростей в НСК 
	double** orbit1 = intergrate(JD, STEP, 6, vec); 
	
	/*
	double** res = orbit1;
	for(int i = 0; i < 1000; i++){
		cout << "time: " << res[i][0] << " x: " << res[i][1] << " y: " << res[i][2] << " z: " << res[i][3] << endl;
		cout << "Vx: " << res[i][3] << " Vy: " << res[i][4] << " Vz: " << res[i][5] << endl;
		cout << endl;	
	}
	*/
	

	/*
	По программе GRACE между спутниками соблюдалось расстояние ~220 +- 50 км
	После построение первой орбиты заметили, что через шаг расстояние от начальной точки примерно такое
	Поэтому за начальную точку второй орбиты возьмём координаты после этого шага
	*/
	vec[0] = 1248.77, vec[1] = -6535.31, vec[2] = 3.71;
	vec[3] = 3.7, vec[4] = 7.49, vec[5] = 1.38;
	
	double **orbit2 = intergrate(JD, STEP, 6, vec);
	
	for(int i = 0; i < 86400.0 / STEP; i++){
		delete[] orbit1[i];
		delete[] orbit2[i];
	}
	delete[] vec;
}

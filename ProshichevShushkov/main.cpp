#include "main.h"

using namespace std;

int main()
{	
	double JD_start = JD;
	double *vec = new double[6];
	int cnt =  86400.0 / STEP;
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
	for(int i = 0; i < cnt; i++){
		cout << "time: " << orbit1[i][0] << " x: " << orbit1[i][1] << " y: " << orbit1[i][2] << " z: " << orbit1[i][3] << endl;
		cout << "Vx: " << orbit1[i][4] << " Vy: " << orbit1[i][5] << " Vz: " << orbit1[i][6] << endl;
		cout << endl;	
	}
	*/
	
	fstream orbit("orbit.txt", ios::out);
	for(int i = 0; i < cnt; i++){
		orbit << orbit1[i][1] << " " << orbit1[i][2] << " " << orbit1[i][3] << endl;
	}
	orbit.close();


	/*
	По программе GRACE между спутниками соблюдалось расстояние ~220 +- 50 км
	После построение первой орбиты заметили, что через шаг расстояние от начальной точки примерно такое
	Поэтому за начальную точку второй орбиты возьмём координаты после этого шага
	*/


	vec[0] = 1472.62, vec[1] = -6718.5, vec[2] = -0.148523;
	vec[3] = 7.43608, vec[4] = 1.63027, vec[5] = 0.000242242;
	
	double **orbit2 = intergrate(JD, STEP, 6, vec);
	
	/*
	fstream file("orbit.txt", ios::out);
	for(int i = 0; i < cnt; i++){
		file << orbit2[i][1] << " " << orbit2[i][2] << " " << orbit2[i][3] << endl;
	}
	file.close();
	*/

	double* res = OrbitDistance(orbit1, orbit2, cnt);
	
	/*
	for(int i = 0; i < cnt; i++){
		cout << res[i][0] << " " << res[i][1] << endl;
	}
	*/

	fstream dist("distance.txt", ios::out);
	for(int i = 0; i < cnt; i++){
		dist << res[2 * i] - JD << " " << res[2 * i + 1] << endl;
	}
	dist.close();


	for(int i = 0; i < cnt; i++){
		delete[] orbit1[i];
		delete[] orbit2[i];
	}
	delete[] vec;			
	delete[] res;
}

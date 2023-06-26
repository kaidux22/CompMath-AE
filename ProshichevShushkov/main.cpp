#include "main.h"

using namespace std;

int main(){
	double JD_start = JD;
	double *vec = new double[12];
<<<<<<< HEAD
	int cnt =  86400.0 / STEP;

=======
	int cnt =  GENERAL_TIME / STEP;
	
>>>>>>> d9a997fa8a9e59c18ef81a71abaa4a30aea08804

	//начальное положение в НСК первого спутника
	vec[0] = 1248.77, vec[1] = -117.887, vec[2] = -6762.66;

	//начальное положение в НСК второго спутника
	vec[3] = 1472.62, vec[4] = -117.105, vec[5] = -6717.48;

	// начальная скорость в НСК первого спутника
	vec[6] = 7.48616, vec[7] = 0.0238816, vec[8] = 1.38195;

	//начальная скорость в НСК второго спутника
	vec[9] = 7.43608, vec[10] = 0.0282099, vec[11] = 1.63003;

	double angle = ANGLE * M_PI / 180.0;

	double rotateMatrix[3][3] = {{1.0, 0.0, 0.0}, {0.0, cos(angle), -sin(angle)}, {0, sin(angle), cos(angle)}};

	for(int i = 0; i < 4; i++){
		changeCoords(rotateMatrix, vec, 3 * i);
	}

	for(int i = 0; i < 12; i++)
		cout << vec[i] << " ";
	cout << endl;

	double** orbits = Integrate(JD, STEP, 12, vec); 

	for(int i = 0; i < 12; i++)
		cout << orbits[0][i + 1] << " ";
	cout << endl;
	
	fstream orbit("orbit.txt", ios::out);
	for(int i = 0; i < cnt; i++){
		orbit << orbits[i][1] << " " << orbits[i][2] << " " << orbits[i][3] << endl;
	}
	orbit.close();

	/*
	По программе GRACE между спутниками соблюдалось расстояние ~220 +- 50 км
	После построение первой орбиты заметили, что через шаг расстояние от начальной точки примерно такое
	Поэтому за начальную точку второй орбиты возьмём координаты после этого шага
	*/

	double* res = OrbitDistance(orbits, cnt);

	
	for(int i = 0; i < cnt; i++){
		res[2 * i + 1] *= (1 + 0 * (rand() % (int)2e5 - 1e5) / 1e8);
	}

	fstream dist("distance.txt", ios::out);
	for(int i = 0; i < cnt; i++){
		dist << res[2 * i] - JD << " " << res[2 * i + 1] << endl;
	}
	dist.close();

	LeastSquare* solve = new LeastSquare(res, cnt);
	solve->Iteration(5);
	delete solve;

	for(int i = 0; i < cnt; i++){
		delete[] orbits[i];
	}
	delete[] vec;			
}

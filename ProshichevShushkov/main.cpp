#include "main.h"

using namespace std;

int main(){
	double JD_start = JD;
	double *vec = new double[12];
	int cnt =  86400.0 / STEP;


	//начальное положение в НСК первого спутника
	vec[0] = 1248.77, vec[1] = -117.887, vec[2] = -6762.66;

	//начальное положение в НСК второго спутника
	vec[3] = 1472.62, vec[4] = -117.105, vec[5] = -6717.48;

	// начальная скорость в НСК первого спутника
	vec[6] = 7.48616, vec[7] = 0.0238816, vec[8] = 1.38195;

	//начальная скорость в НСК второго спутника
	vec[9] = 7.43608, vec[10] = 0.0282099, vec[11] = 1.63003;

	double** orbits = Integrate(JD, STEP, 12, vec); 
	

	/*
	for(int i = 0; i < cnt; i++){
		cout << "time: " << orbit1[i][0] << " x: " << orbit1[i][1] << " y: " << orbit1[i][2] << " z: " << orbit1[i][3] << endl;
		cout << "Vx: " << orbit1[i][4] << " Vy: " << orbit1[i][5] << " Vz: " << orbit1[i][6] << endl;
		cout << endl;	
	}
	*/
	
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
	
	/*
	fstream file("orbit.txt", ios::out);
	for(int i = 0; i < cnt; i++){
		file << orbit2[i][1] << " " << orbit2[i][2] << " " << orbit2[i][3] << endl;
	}
	file.close();
	*/

	double* res = OrbitDistance(orbits, cnt);
	
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

	LeastSquare* solve = new LeastSquare(res, cnt);
	solve->Iteration(5);
	delete solve;

	for(int i = 0; i < cnt; i++){
		delete[] orbits[i];
	}
	delete[] vec;			
}

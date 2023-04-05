#include "main.h"

#define JD 2451545.0
#define STEP 60.0

using namespace std;

//f(x(t)) = -k*x(t)

//-4661,36 -3953,12 3154,59 
//-4742,49 -4021,92 2939,37

int main()
{	double JD_start = JD;
	double *vec = new double[3];
	double rotateMatrix[3][3];
	vec[0] = -4661.36, vec[1] = -3953.12, vec[2] = 3154.59;
	iauC2t06a(JD_start + (37.0 + 32.184) / 86400.0, 0, JD_start, 0, 0, 0, rotateMatrix);

	Transposition(rotateMatrix);
	changeCoords(rotateMatrix, vec);
	double** res = intergrate(JD, STEP, 3, vec); 
	
	for(int i = 0; i < 1000; i++){
		cout << "time: " << res[i][0] << " x: " << res[i][1] << " y: " << res[i][2] << " z: " << res[i][3] << endl;
	}
	delete[] vec;
}

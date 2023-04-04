// CompMath.cpp: определяет точку входа для приложения.
//

#include "main.h"

#define JD 86400.0

using namespace std;

//f(x(t)) = -k*x(t)

//-4661,36 -3953,12 3154,59 
//-4742,49 -4021,92 2939,37

int main()
{	double UTC_start = JD * 24.0 * 60.0 * 60.0;
	double *vec = new double[3];
	double rotateMatrix[3][3];
	vec[0] = -4661.36, vec[1] = -3953.12, vec[2] = 3154.59;
	iauC2t06a(UTC_start + 37 + 32.184, 0, UTC_start, 0, 0, 0, rotateMatrix);
	
	Transposition(rotateMatrix);
	changeCoords(rotateMatrix, vec);
	
	double *res = GradV(vec, UTC_start);
	cout << res[0] << " " << res[1] << " " << res[2] << endl;
}

// CompMath.cpp: определяет точку входа для приложения.
//

#include "main.h"

using namespace std;

//f(x(t)) = -k*x(t)

//-4661,36 -3953,12 3154,59 
//-4742,49 -4021,92 2939,37

int main()
{
	double rotateMatrix[3][3];
	iauC2t06a(2451545.0 * 24.0 * 60.0 * 60.0 + 37.0 + 32.184, 0, 2451545.0 * 24.0 * 60.0 * 60.0, 0, 0, 0, rotateMatrix);

	for(int i = 0; i < 3; i++){
		for(int j = 0; j < 3; j++){
			cout << rotateMatrix[i][j] << " ";
		}
		cout << endl;
	}
}

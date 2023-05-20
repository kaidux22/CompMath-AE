
#include <utility>
#include <vector>

using namespace std;

/*
Функция транспонирует матрицу 3х3
*/
void Transposition(double(*matrix)[3]);

/*
Cмена координат
*/
void changeCoords(double(*rotateMatrix)[3], double* vec, int idx);
void changeCoordsRight(double(*rotateMatrix)[3], double* vec, int idx);

double** multiplication_AtA(vector<vector<double>>& A);
double* multiplication_Atr(vector<vector<double>> &A, vector<double>& r);
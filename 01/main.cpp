#include "matrix.h"
#include <fstream>
#include <iostream>
#include <random>

//using namespace std;

void fill_matrixs(Matrix &X, Matrix& Y)
{
  std::ifstream data("ACR.txt");
  for (int i = 0; i < X.GetNumRows(); i++)
  {
    for (int j = 0; j < X.GetNumColumns(); j++)
      data >> X[i][j];
    data >> Y[i][0];
  }
}
vector<double> int_coef(vector<double> yt, const vector<double>& x, const Matrix& O);

int main()
{
  Matrix X(100, 7);
  Matrix Y(100, 1);
  fill_matrixs(X, Y);
  Matrix result = MatrixE(Transp(X)*X) * Transp(X) * Y;
  Matrix P(7, 7);
  for (int i = 0; i < 7; i++)
    P[i][i] = 1'000'000;
  Matrix O(7, 1);
  for (int i = 0; i < 100; i++) {
    P = P + (-1.0) * (P * (X[i] * X[i]) * P) / (1 + (X[i] * P * X[i])[0][0]);
    O = O + (P * X[i]) * int_coef(Y[i], X[i], O);
  }
  std::cout << /*"P:\n" << P <<*/ "Omega:\n" << O << "Etalon:\n" << result;
}

vector<double> int_coef(vector<double> yt, const vector<double>& x, const Matrix& O) {
  Matrix xO = x*O;
  for (int i = 0; i < yt.size(); i++)
    yt[i] -= xO[0][i];
  return yt;
}
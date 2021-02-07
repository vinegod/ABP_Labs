#include "matrix.h"
#include <fstream>
#include <iostream>
#include <ctime>
//#include <random>

//using namespace std;

void create_data(const int N) {
  const Matrix params({ {0}, {0.22}, {-0.18}, {-0.8}, {1}, {0.5}, {0.5}, {0} });

  Matrix X(N, params.GetNumRows());
  X[0][0] = 1;
  for (int j = 1; j < X.GetNumColumns(); j++)
    X[0][j] = (rand()%50 - 25); //Первые полностью рандомные значения
  Matrix Y(N, 1);
  Y[0][0] = (column_matrix(X[0]) * params)[0][0];

  for (int i = 1; i < N; i++) {
    X[i][0] = 1;
    X[i][1] = Y[i - 1][0];
    int j = 2;
    for (j = 2; j < X.GetNumColumns()/2; j++)
      X[i][j] = X[i - 1][j - 1];
    X[i][j] = (rand()%5 - 2.5) * 1.0 / j;
    for (j++; j < X.GetNumColumns(); j++)
      X[i][j] = X[i - 1][j - 1];
    Y[i][0] = (column_matrix(X[i]) * params)[0][0];
  }
  std::ofstream data("generated_data.txt");
  for (int i = 0; i < X.GetNumRows(); i++)
  {
    for (int j = 0; j < X.GetNumColumns(); j++)
      data << X[i][j] <<'\t';
    data << Y[i][0] << '\n';
  }
  data.close();
}

void fill_matrixs(Matrix &X, Matrix& Y)
{
  std::ifstream data("generated_data.txt");
  for (int i = 0; i < X.GetNumRows(); i++)
  {
    for (int j = 0; j < X.GetNumColumns(); j++)
      data >> X[i][j];
    data >> Y[i][0];
  }
  data.close();
}
vector<double> int_coef(vector<double> yt, const vector<double>& x, const Matrix& O);

int main()
{
  srand(time(NULL));
  const int N = 100;
  const int params_count = 8;
  create_data(N);
  Matrix X(N, params_count);
  Matrix Y(N, 1);
  fill_matrixs(X, Y);
  Matrix result = MatrixE(Transp(X)*X) * Transp(X) * Y;
  std::cout << result;
  Matrix P(params_count, params_count);
  for (int i = 0; i < params_count; i++)
    P[i][i] = 1'000'000;
  Matrix O(params_count, 1);
  for (int i = 0; i < N; i++) {
    P = P + (-1.0) * (P * (X[i] * X[i]) * P) / (1 + (X[i] * P * X[i])[0][0]);
    O = O + (P * X[i]) * int_coef(Y[i], X[i], O);
  }
  std::cout << "Omega:\n" << O << "Etalon:\n" << result;
}

vector<double> int_coef(vector<double> yt, const vector<double>& x, const Matrix& O) {
  Matrix xO = x*O;
  for (int i = 0; i < yt.size(); i++)
    yt[i] -= xO[0][i];
  return yt;
}
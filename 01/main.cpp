#include "matrix.h"
#include <fstream>
#include <iostream>
#include <ctime>
//#include <random>

//using namespace std;

void create_data(const int N);
void fill_matrixs(Matrix &X, Matrix& Y);
Matrix count_omega(Matrix&, Matrix&, int params_count, int betta);

int main()
{
  srand(time(NULL));
  const int N = 1000;
  const int params_count = 8;
  //create_data(N);
  Matrix X(N, params_count);
  Matrix Y(N, 1);
  fill_matrixs(X, Y);
  Matrix result = MatrixE(Transp(X)*X) * Transp(X) * Y;
  Matrix O = count_omega(X, Y, params_count, 1000);
  std::cout << "Omega:\n" << O << "Etalon:\n" << result;
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

Matrix count_omega(Matrix& X, Matrix& Y, int params_count, int betta){
  Matrix P(params_count, params_count);
  for (int i = 0; i < params_count; i++)
    P[i][i] = betta;
  Matrix O(params_count, 1);
  for (int i = 0; i < X.GetNumRows(); i++) {
    P = P + (-1.0) * (P * row_matrix(X[i]) * column_matrix(X[i]) * P) /
      (1 + (column_matrix(X[i]) * P * row_matrix(X[i]))[0][0]);
    O = O + P * row_matrix(X[i]) * (row_matrix(Y[i]) + (-1.0) * column_matrix(X[i]) * O);
  }
  return O;
}
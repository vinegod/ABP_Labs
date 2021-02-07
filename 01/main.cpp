#include "matrix.h"
#include <cmath>
#include <ctime>
#include <fstream>
#include <iostream>

//#include <random>

// using namespace std;

void create_data(const int N);
void fill_matrixs(Matrix &X, Matrix &Y);
Matrix count_omega(Matrix &X, Matrix &Y, int params_count, int betta);
std::pair<double, double> Akaiki(Matrix &, Matrix &, Matrix &);

int main() {
  srand(time(NULL));
  const int N = 1000;
  const int params_count = 8;
  create_data(N);
  Matrix X(N, params_count);
  Matrix Y(N, 1);
  fill_matrixs(X, Y);
  Matrix MHK = MatrixE(Transp(X) * X) * Transp(X) * Y;
  Matrix PMHK = count_omega(X, Y, params_count, 10'000'000);
  std::cout << "MHK algorithm:\n" << MHK << "PMHK algorithm:\n" << PMHK;

  auto Akaiki_MHK = Akaiki(X, Y, MHK);
  std::cout << "Sum E^2 for MHK: " << Akaiki_MHK.first << '\n'
            << "Akaiki for MHK: " << Akaiki_MHK.second << '\n';
  auto Akaiki_PMHK = Akaiki(X, Y, PMHK);
  std::cout << "Sum E^2 for PMHK: " << Akaiki_PMHK.first << '\n'
            << "Akaiki for PMHK: " << Akaiki_PMHK.second << '\n';
}

void fill_matrixs(Matrix &X, Matrix &Y) {
  std::ifstream data("generated_data.txt");
  for (int i = 0; i < X.GetNumRows(); i++) {
    for (int j = 0; j < X.GetNumColumns(); j++)
      data >> X[i][j];
    data >> Y[i][0];
  }
  data.close();
}

void create_data(const int N) {
  const Matrix params({{0},
                       {0.22},
                       {-0.18},
                       {0.08},
                       {1},
                       {0.5},
                       {0.5},
                       {0}}); //Параметры, которые были заданы
  const uint16_t count_A = 3;
  //const uint16_t count_B = 3;

  Matrix X(N, params.GetNumRows());
  X[0][0] = 1;
  for (int j = 1; j < X.GetNumColumns(); j++)
    X[0][j] = (rand() % 5 - 2.5); //Первые полностью рандомные значения
  Matrix Y(N, 1);
  Y[0][0] = (column_matrix(X[0]) * params)[0][0];

  for (int i = 1; i < N; i++) {
    X[i][0] = 1;
    X[i][1] = Y[i - 1][0];
    int j;
    for (j = 2; j < 1 + count_A; j++)
      X[i][j] = X[i - 1][j - 1];
    X[i][j] = (rand() % 5 - 2.5) * 1.0 / j;
    for (++j; j < X.GetNumColumns(); j++)
      X[i][j] = X[i - 1][j - 1];
    Y[i][0] = (column_matrix(X[i]) * params)[0][0];
  }
  std::ofstream data("generated_data.txt");
  for (int i = 0; i < X.GetNumRows(); i++) {
    for (int j = 0; j < X.GetNumColumns(); j++)
      data << X[i][j] << '\t';
    data << Y[i][0] << '\n';
  }
  data.close();
}

Matrix count_omega(Matrix &X, Matrix &Y, int params_count, int betta) {
  Matrix P(params_count, params_count);
  for (int i = 0; i < params_count; i++)
    P[i][i] = betta;
  Matrix O(params_count, 1);
  for (int i = 0; i < X.GetNumRows(); i++) {
    P = P + (-1.0) * (P * row_matrix(X[i]) * column_matrix(X[i]) * P) /
                (1 + (column_matrix(X[i]) * P * row_matrix(X[i]))[0][0]);
    O = O + P * row_matrix(X[i]) *
                (row_matrix(Y[i]) + (-1.0) * column_matrix(X[i]) * O);
  }
  return O;
}

std::pair<double, double> Akaiki(Matrix &X, Matrix &Y, Matrix &params) {
  int N = X.GetNumRows();
  double IKA = 0;
  for (int i = 0; i < N; i++) {
    IKA += pow((column_matrix(X[i]) * params + (-1.0) * row_matrix(Y[i]))[0][0], 2);
  }

  return {IKA, (N * log(IKA) / log(2.7182818284) + 2 * N)};
}
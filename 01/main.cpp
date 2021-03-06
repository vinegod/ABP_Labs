#include "matrix.h"
#include <cmath>
#include <ctime>
#include <fstream>
#include <iostream>
#include <tuple>

void create_data(const int N, uint16_t, uint16_t, bool random_data);
void search_params(int count_A, int count_B, int, bool);

int main() {
  //srand(time(NULL));
  const int N = 100;
  const int betta = 1'000;
  //create_data(N, 3, 3, 0);
  /*for (int i = 1; i < 4; i++)
    for (int j = 1; j < 4; j++) {
      std::cout << "For APKC(" << i << ", " << j << "):\n";
      search_params(i, j, betta);
      std::cout << '\n';
    }*/
  search_params(3, 3, N, 1);
}

int create_and_fill_from_file(int count_A, int count_B, Matrix &, Matrix &Y);
void fill_X(Matrix &X, const Matrix &Eps, const Matrix &Y, int, int);
Matrix change_Y_size(int, const Matrix &, bool);
Matrix count_omega(const Matrix &X, const Matrix &Y, int params_count,
                   int betta);
std::tuple<double, double, double> Akaiki(const Matrix &, const Matrix &,
                                          const Matrix &);

Matrix recreate_matrix(const Matrix &old, int size) {
  Matrix new_m(size, old.GetNumColumns());
  for (int i = 0; i < size; i++)
    new_m[i] = old[i];
  return new_m;
}
void vector_diff(const Matrix &X, const Matrix &Y) {

  const Matrix params =
      Matrix({{0}, {0.22}, {-0.18}, {0.08}, {1}, {0.5}, {0.5}, {0}});
  std::ofstream fk1("stupid_data1.txt");
  std::ofstream fk2("stupid_data2.txt");
  //fk << "Iter" << '\t' << "Norm_MHK" << '\t' << "Norm_PMNK" << '\n';
  for (int i = 2; i < X.GetNumRows(); i++) {
    Matrix new_X = recreate_matrix(X, i);
    Matrix new_Y = recreate_matrix(Y, i);
    Matrix MHK = MatrixE(Transp(new_X) * new_X) * Transp(new_X) * new_Y;
    Matrix PMHK = count_omega(new_X, new_Y, params.GetNumRows(), 1'000'000);
    fk1 << Transp(MHK);
    fk2 << Transp(PMHK);
  }
  //std::cout << count_omega(X, Y, params.GetNumRows(), 1'000'000);
}

void search_params(int count_A, int count_B, int betta, bool vectors_diff) {
  const int params_count = count_A + count_B + 2;
  Matrix Eps, Y;
  int N = create_and_fill_from_file(count_A, count_B, Eps, Y);
  Matrix X(N, params_count);
  fill_X(X, Eps, Y, count_A, count_B);
  Y = change_Y_size(N, Y, count_A > count_B);
  
  if (vectors_diff)
    vector_diff(X, Y);
  else {
    Matrix MHK = MatrixE(Transp(X) * X) * Transp(X) * Y;
    Matrix PMHK = count_omega(X, Y, params_count, betta);

    std::cout << "MHK algorithm parameters:\n"
              << MHK << "PMHK algorithm parameters:\n"
              << PMHK;

    double MHK_S, MHK_R, MHK_E;
    std::tie(MHK_S, MHK_R, MHK_E) = Akaiki(X, Y, MHK);
    std::cout << "\nSum E^2 for MHK: " << MHK_S << '\n'
              << "R^2 for MHK: " << MHK_R << '\n'
              << "Akaiki for MHK: " << MHK_E << '\n';
    /*std::cout << MHK_S << '\t'
              << MHK_R << '\t'
              << MHK_E << '\n';*/

    double PMHK_S, PMHK_R, PMHK_E;
    std::tie(PMHK_S, PMHK_R, PMHK_E) = Akaiki(X, Y, PMHK);
    std::cout << "Sum E^2 for MHK: " << PMHK_S << '\n'
              << "R^2 for MHK: " << PMHK_R << '\n'
              << "Akaiki for MHK: " << PMHK_E << '\n';
    /*std::cout << PMHK_S << '\t'
              << PMHK_R << '\t'
              << PMHK_E << '\n';*/
  }
}

int create_and_fill_from_file(int count_A, int count_B, Matrix &Eps,
                              Matrix &Y) {
  std::ifstream Eps_data("Eps_data.txt");
  std::ifstream Y_data("Y_data.txt");

  vector<double> y;
  double temp;
  while (Y_data >> temp) {
    y.emplace_back(temp);
  }
  vector<double> eps;
  while (Eps_data >> temp) {
    eps.emplace_back(temp);
  }
  Eps = row_matrix(eps);
  Y = row_matrix(y);

  Eps_data.close();
  Y_data.close();
  return std::min(y.size(), eps.size()) - std::max(count_A, count_B + 1);
}

void fill_X(Matrix &X, const Matrix &Eps, const Matrix &Y, int count_A,
            int count_B) {
  int n = std::max(count_A, count_B + 1);
  X[0][0] = 1;
  if (count_A > count_B)
    n++;
  for (int j = 1; j < 1 + count_A; j++)
    X[0][j] = Y[n - j - 1][0];

  for (int j = 1 + count_A; j < X.GetNumColumns(); j++)
    X[0][j] = Eps[n + count_A - j][0];

  for (int i = 1; i < X.GetNumRows(); i++) {
    X[i][0] = 1;
    X[i][1] = Y[i + n - 2][0];
    for (int j = 2; j < 1 + count_A; j++)
      X[i][j] = X[i - 1][j - 1];
    X[i][1 + count_A] = Eps[i + n - 1][0];
    for (int j = 2 + count_A; j < X.GetNumColumns(); j++)
      X[i][j] = X[i - 1][j - 1];
  }
  // std::ofstream data_X("X.txt");
  // data_X << X;
}

Matrix change_Y_size(int N, const Matrix &Y, bool is_less) {
  Matrix new_Y(N, 1);
  for (int i = 0; i < N; i++) {
    if (is_less)
      new_Y[i][0] = Y[Y.GetNumRows() - N + i][0];
    else
      new_Y[i][0] = Y[Y.GetNumRows() - N + i - 1][0];
  }
  // std::ofstream data_Y("Y.txt");
  // data_Y << new_Y;
  return new_Y;
}

Matrix count_omega(const Matrix &X, const Matrix &Y, int params_count,
                   int betta) {
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

std::tuple<double, double, double> Akaiki(const Matrix &X, const Matrix &Y,
                                          const Matrix &params) {
  int N = X.GetNumRows();
  double SS_res = 0, SS_tot = 0;
  double y_ = 0;
  for (int i = 0; i < N; i++)
    y_ += Y[i][0] / N;
  for (int i = 0; i < N; i++) {
    SS_res += pow(
        (column_matrix(X[i]) * params + (-1.0) * row_matrix(Y[i]))[0][0], 2);
    SS_tot += pow((Y[i][0] - y_), 2);
  }

  return std::make_tuple(
      SS_res, 1 - SS_res / SS_tot,
      (N * log(SS_res) / log(2.7182818284) + 2 * params.GetNumRows()));
}

void create_data(const int N, uint16_t count_A, uint16_t count_B,
                 bool random_data) {
  Matrix params;
  if (random_data) {
    vector<double> parameters;
    for (uint16_t i = 0; i < count_A + 1; i++)
      parameters.push_back((rand() % 1000 - 500) * 1.0 / 1000);
    parameters.push_back(1);
    for (uint16_t i = 0; i < count_B; i++)
      parameters.push_back((rand() % 1000 - 500) * 1.0 / 1000);
    params = row_matrix(parameters);
    std::cout << "APKC(" << count_A << ", " << count_B << "):\n"
              << params << '\n';
  } else {
    params = Matrix({{0},
                     {0.22},
                     {-0.18},
                     {0.08},
                     {1},
                     {0.5},
                     {0.5},
                     {0}}); //Параметры, которые были заданы
  }
  Matrix X(N, params.GetNumRows());
  X[0][0] = 1;

  for (int j = 1; j < X.GetNumColumns(); j++) {
    X[0][j] =
        (rand() % 100 - 50) * 1.0 / 100; //Первые полностью рандомные значения
  }
  Matrix Y(N, 1);
  Y[0][0] = (column_matrix(X[0]) * params)[0][0];
  for (int i = 1; i < N; i++) {
    X[i][0] = 1;
    X[i][1] = Y[i - 1][0];
    int j;
    for (j = 2; j < 1 + count_A; j++)
      X[i][j] = X[i - 1][j - 1];
    X[i][j] = (rand() % 100 - 50) * 1.0 / 100;
    for (++j; j < X.GetNumColumns(); j++)
      X[i][j] = X[i - 1][j - 1];
    Y[i][0] = (column_matrix(X[i]) * params)[0][0];
  }
  std::ofstream Eps_data("Eps_data.txt");
  std::ofstream X_data("X_data.txt");
  std::ofstream Y_data("Y_data.txt");
  for (int i = 0; i < X.GetNumRows(); i++) {
    for (int j = 0; j < 1 + count_A; j++)
      X_data << X[i][j] << '\t';
    for (int j = 1 + count_A; j < X.GetNumColumns(); j++)
      X_data << X[i][j] << '\t';
    Eps_data << X[i][1 + count_A] << '\n';
    Y_data << Y[i][0] << '\n';
    X_data << Y[i][0] << '\n';
  }
  X_data.close();
  Eps_data.close();
  Y_data.close();
}
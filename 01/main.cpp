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

int main()
{
  Matrix X(100, 7);
  Matrix Y(100, 1);
  fill_matrixs(X, Y);
  Matrix XT = Transp(X);
  Matrix XTX = XT*X;
  Matrix result = MatrixE(XTX) * XT * Y;
  std::cout << result;

}

void func()
{
  int p = 3, q = 2;
  std::default_random_engine generator;
  std::normal_distribution<double> distribution(5.0, 2.0);
  double mass[10] = {};
  for (int i = 0; i < 1000; ++i)
  {
    double number = distribution(generator);
    if ((number >= 0.0) && (number < 10.0))
      ++mass[int(number)];
  }
  vector<double> a({mass[0], mass[1], mass[2]});
  for (auto &v : mass)
    std::cout << v << ' ';
}
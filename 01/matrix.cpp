#include "matrix.h"
#include <cmath>

Matrix::Matrix() {
  num_rows = 0;
  num_columns = 0;
}

Matrix::Matrix(vector<vector<double>> &&elements_) {
  elements = elements_;
  num_rows = elements.size();
  num_columns = elements[0].size();
}

Matrix::Matrix(int num_rows, int num_columns) { Reset(num_rows, num_columns); }

void Matrix::Reset(int num_rows_, int num_columns_) {

  num_rows = num_rows_;
  num_columns = num_columns_;
  elements.assign(num_rows, vector<double>(num_columns));
}

int Matrix::GetNumRows() const { return num_rows; }

int Matrix::GetNumColumns() const { return num_columns; }

bool operator==(const Matrix &one, const Matrix &two) {
  if (one.GetNumRows() != two.GetNumRows()) {
    return false;
  }

  if (one.GetNumColumns() != two.GetNumColumns()) {
    return false;
  }

  for (int row = 0; row < one.GetNumRows(); row++) {
    for (int column = 0; column < one.GetNumColumns(); column++) {
      if (one[row][column] != two[row][column]) {
        return false;
      }
    }
  }

  return true;
}

Matrix operator+(const Matrix &one, const Matrix &two) {
  if (one.GetNumRows() != two.GetNumRows()) {
    throw std::invalid_argument("Mismatched number of rows");
  }

  if (one.GetNumColumns() != two.GetNumColumns()) {
    throw std::invalid_argument("Mismatched number of columns");
  }

  Matrix result(one.GetNumRows(), one.GetNumColumns());
  for (int row = 0; row < result.GetNumRows(); row++) {
    for (int column = 0; column < result.GetNumColumns(); column++) {
      result[row][column] = one[row][column] + two[row][column];
    }
  }

  return result;
}

Matrix operator*(const Matrix &one, const Matrix &two) {
  if (one.GetNumColumns() != two.GetNumRows())
    throw std::invalid_argument("Mismatched number of rows and columns");
  Matrix result(one.GetNumRows(), two.GetNumColumns());
  for (int i = 0; i < one.GetNumRows(); i++)
    for (int k = 0; k < one.GetNumColumns(); k++)
      for (int j = 0; j < two.GetNumColumns(); j++)
        result[i][j] += one[i][k] * two[k][j];
  return result;
}
Matrix operator/(const Matrix &matrix, double x) {
  Matrix result(matrix.GetNumRows(), matrix.GetNumColumns());
  for (int i = 0; i < matrix.GetNumRows(); i++)
    for (int j = 0; j < matrix.GetNumColumns(); j++)
      result[i][j] = matrix[i][j] / x;
  return result;
}

Matrix Transp(const Matrix &matrix) {
  Matrix result(matrix.GetNumColumns(), matrix.GetNumRows());
  for (int i = 0; i < matrix.GetNumColumns(); i++)
    for (int j = 0; j < matrix.GetNumRows(); j++)
      result[i][j] = matrix[j][i];
  return result;
}

std::istream &operator>>(std::istream &in, Matrix &matrix) {
  int num_rows, num_columns;
  in >> num_rows >> num_columns;

  matrix.Reset(num_rows, num_columns);
  for (int row = 0; row < num_rows; row++) {
    for (int column = 0; column < num_columns; column++) {
      in >> matrix[row][column];
    }
  }

  return in;
}

const std::vector<double> &Matrix::operator[](size_t i) const {
  return elements[i];
}
std::vector<double> &Matrix::operator[](size_t i) { return elements[i]; }

std::ostream &operator<<(std::ostream &out, const Matrix &matrix) {
  for (int row = 0; row < matrix.GetNumRows(); row++) {
    for (int column = 0; column < matrix.GetNumColumns(); column++) {
      if (column > 0) {
        out << " ";
      }
      out << matrix[row][column];
    }
    out << std::endl;
  }
  return out;
}

Matrix operator*(double one, const Matrix &two) {
  Matrix result(two.GetNumRows(), two.GetNumColumns());
  for (uint16_t i = 0; i < two.GetNumRows(); i++)
    for (int j = 0; j < two.GetNumColumns(); j++)
      result[i][j] = one * two[i][j];

  return result;
}

void Gauss(Matrix &matrix) {
  for (int i = 0; i < matrix.GetNumRows(); i++) {
    int row = i;
    int mx = matrix[i][i];
    for (int k = i + 1; k < matrix.GetNumRows(); k++) {
      if (fabs(matrix[k][i]) > mx) {
        row = k;
        mx = fabs(matrix[k][i]);
      }
    }
    if (row != i)
      std::swap(matrix[row], matrix[i]);
    double first = matrix[i][i];
    for (int k = i; k < matrix.GetNumColumns(); k++)
      matrix[i][k] /= first;
    for (int j = i + 1; j < matrix.GetNumRows(); j++) {
      double e = matrix[j][i] / matrix[i][i];
      for (int k = i; k < matrix.GetNumColumns(); k++)
        matrix[j][k] -= e * matrix[i][k];
    }
  }
  for (int i = matrix.GetNumRows() - 1; i >= 0; i--)
    for (int j = i - 1; j >= 0; j--) {
      double e = matrix[j][i] / matrix[i][i];
      for (int k = matrix.GetNumColumns() - 1; k >= 0; k--)
        matrix[j][k] -= e * matrix[i][k];
    }
}

Matrix MatrixE(const Matrix &matrix) {
  Matrix temp(matrix.GetNumRows(), 2 * matrix.GetNumColumns());
  for (int i = 0; i < matrix.GetNumRows(); i++)
    for (int j = 0; j < matrix.GetNumColumns(); j++)
      temp[i][j] = matrix[i][j];
  for (int i = matrix.GetNumColumns(); i < 2 * matrix.GetNumColumns(); i++)
    temp[i - matrix.GetNumColumns()][i] = 1;
  Gauss(temp);
  Matrix result(matrix.GetNumRows(), matrix.GetNumColumns());
  for (int i = 0; i < matrix.GetNumRows(); i++)
    for (int j = 0; j < matrix.GetNumColumns(); j++)
      result[i][j] = temp[i][matrix.GetNumColumns() + j];
  return result;
}

Matrix row_matrix(const vector<double> &v) {
  Matrix result(v.size(), 1);
  for (uint16_t i = 0; i < v.size(); i++)
    result[i][0] = v[i];
  return result;
}
Matrix column_matrix(const vector<double> &v) {
  Matrix result(1, v.size());
  for (uint16_t i = 0; i < v.size(); i++)
    result[0][i] = v[i];
  return result;
}
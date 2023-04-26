#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>

using namespace std;

class Matrix {
private:
    int rows, columns, size = 0;
    vector<vector<double>> data;

public:
    Matrix(int row, int column) {
        this->rows = row;
        this->columns = column;
        this->size = row;
        this->data.resize(row, vector<double>(column));
    }

    int get_rows() const {
        return rows;
    }

    int get_columns() const {
        return columns;
    }

    virtual int get_size() const {
        return size;
    }

    virtual void set_size(int n) {
        size = n;
    }

    double& operator()(int i, int j) {
        return data[i][j];
    }

    const double& operator()(int i, int j) const {
        return data[i][j];
    }

    Matrix operator+(const Matrix& secMatrix) const {
        if (rows != secMatrix.rows || columns != secMatrix.columns) {
            cout << "Error: the dimensional problem occurred\n";
            Matrix result(0, 0);
            return result;
        }
        Matrix result(rows, columns);
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < columns; j++) {
                result(i, j) = (*this)(i, j) + secMatrix(i, j);
            }
        }
        return result;
    }

    Matrix operator-(const Matrix& secMatrix) const {
        if (rows != secMatrix.rows || columns != secMatrix.columns) {
            cout << "Error: the dimensional problem occurred\n";
            Matrix result(0, 0);
            return result;
        }
        Matrix result(rows, columns);
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < columns; j++) {
                result(i, j) = (*this)(i, j) - secMatrix(i, j);
            }
        }
        return result;
    }

    Matrix operator=(const Matrix* secMatrix) const {
        (*this) = secMatrix;
        return *secMatrix;
    }

    Matrix operator*(const Matrix& secMatrix) const {
        if (columns != secMatrix.rows) {
            cout << "Error: the dimensional problem occurred\n";
            Matrix result(0, 0);
            return result;
        }
        Matrix result(rows, secMatrix.columns);
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < secMatrix.columns; j++) {
                double sum = 0;
                for (int k = 0; k < columns; k++) {
                    sum += (*this)(i, k) * secMatrix(k, j);
                }
                result(i, j) = sum;
            }
        }
        return result;
    }

    Matrix transpose() const {
        Matrix result(columns, rows);
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < columns; j++) {
                result(j, i) = (*this)(i, j);
            }
        }
        return result;
    }

    double determinant() const {
        if (rows != columns) {
            cout << "Error: the dimensional problem occurred\n";
            return 0;
        }
        if (rows == 1) {
            return (*this)(0, 0);
        }
        double result = 0;
        for (int i = 0; i < rows; i++) {
            Matrix temp(rows - 1, columns - 1);
            for (int j = 1; j < rows; j++) {
                for (int k = 0; k < columns; k++) {
                    if (k < i) {
                        temp(j - 1, k) = (*this)(j, k);
                    }
                    else if (k > i) {
                        temp(j - 1, k - 1) = (*this)(j, k);
                    }
                }
            }
            result += pow(-1, i) * (*this)(0, i) * temp.determinant();
        }
        cout << result << "\n";
        return result;
    }

    Matrix inverse() const {
        if (rows != columns) {
            cout << "Error: the dimensional problem occurred\n";
            Matrix result(0, 0);
            return result;
        }
        Matrix augmented(rows, 2 * columns);
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < columns; j++) {
                augmented(i, j) = (*this)(i, j);
            }
            for (int j = columns; j < 2 * columns; j++) {
                augmented(i, j) = (i == j - columns) ? 1.0 : 0.0;
            }
        }
        for (int j = 0; j < columns; j++) {
            int max_row = j;
            double max_val = augmented(j, j);
            for (int i = j + 1; i < rows; i++) {
                if (abs(augmented(i, j)) > abs(max_val)) {
                    max_row = i;
                    max_val = augmented(i, j);
                }
            }
            if (max_val == 0) {
                cout << "Error: the matrix is singular\n";
                Matrix result(0, 0);
                return result;
            }
            if (max_row != j) {
                for (int k = j; k < 2 * columns; k++) {
                    swap(augmented(j, k), augmented(max_row, k));
                }
            }
            double pivot = augmented(j, j);
            for (int k = j; k < 2 * columns; k++) {
                augmented(j, k) /= pivot;
            }
            for (int i = 0; i < rows; i++) {
                if (i != j) {
                    double factor = augmented(i, j);
                    for (int k = j; k < 2 * columns; k++) {
                        augmented(i, k) -= factor * augmented(j, k);
                    }
                }
            }
        }
        Matrix result(rows, columns);
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < columns; j++) {
                result(i, j) = augmented(i, j + columns);
            }
        }
        return result;
    }


};

ostream& operator<<(ostream& outputStream, const Matrix& matrix) {
    for (int i = 0; i < matrix.get_rows(); i++) {
        for (int j = 0; j < matrix.get_columns(); j++) {
            if (j == matrix.get_columns() - 1){
            outputStream << fixed << setprecision(4) << matrix(i, j);
                
            } else {
                outputStream << fixed << setprecision(4) << matrix(i, j) << " ";
            }
        }
        outputStream << endl;
    }
    return outputStream;
}

istream& operator>>(istream& inputStream, Matrix& matrix) {
    for (int i = 0; i < matrix.get_rows(); i++) {
        for (int j = 0; j < matrix.get_columns(); j++) {
            inputStream >> matrix(i, j);
        }
    }
    return inputStream;
}


int main() {
    int size;
    cin >> size;

    int* t = new int[size];
    int* b = new int[size];

    for (int i = 0; i < size; i++) {
        cin >> t[i] >> b[i];
    }

    int polinomial_degree;
    cin >> polinomial_degree;

    Matrix A(size, polinomial_degree + 1);
    Matrix B(size, 1);

    for (int i = 0; i < size; i++) {
        for (int j = 0; j < polinomial_degree + 1; j++) {
            A(i, j) = pow(t[i], j);
        }
        B(i, 0) = b[i];
    }

    cout << "A:\n" << A;

    Matrix A_transpose = A.transpose();

    cout << "A_T*A:\n" << A_transpose * A;

    Matrix A_inverse = (A_transpose * A).inverse();

    cout << "(A_T*A)^-1:\n" << A_inverse;

    Matrix A_transpose_B = A_transpose * B;

    cout << "A_T*b:\n" << A_transpose_B;

    Matrix x = A_inverse * A_transpose_B;

    cout << "x~:\n" << x;

    return 0;
}

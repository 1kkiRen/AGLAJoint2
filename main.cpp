#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <stdio.h>

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

    Matrix operator*(const double& scalar) const {
        Matrix result(rows, columns);
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < columns; j++) {
                result(i, j) = (*this)(i, j) * scalar;
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

class IdentityMatrix : public Matrix {
public:
    IdentityMatrix(int n) : Matrix(n, n) {
        for (int i = 0; i < n; i++) {
            (*this)(i, i) = 1;
        }
    }
};

class PermutationMatrix : public Matrix {
public:
    PermutationMatrix(int n, int i, int j) : Matrix(n, n) {
        for (int k = 0; k < n; k++) {
            (*this)(k, k) = 1;
        }
        (*this)(i, i) = 0;
        (*this)(j, j) = 0;
        (*this)(i, j) = 1;
        (*this)(j, i) = 1;
    }

};

class AugmentedMatrix : public Matrix {
public:
    Matrix* Right;
    Matrix* Left;
    AugmentedMatrix(Matrix A, Matrix B) : Matrix(A.get_size(), A.get_size() + B.get_size()) {
        for (int i = 0; i < A.get_size(); i++) {
            for (int j = 0; j < A.get_size(); j++) {
                (*this)(i, j) = A(i, j);
            }
        }
        for (int i = 0; i < B.get_size(); i++) {
            for (int j = 0; j < B.get_size(); j++) {
                (*this)(i, j + A.get_size()) = B(i, j);
            }
        }
        Left = new Matrix(B.get_size(), B.get_size());
        Right = new Matrix(A.get_size(), A.get_size());
        for (int i = 0; i < B.get_size(); i++) {
            for (int j = 0; j < B.get_size(); j++) {
                Left->operator()(i, j) = (*this)(i, j + A.get_size());
            }
        }
        for (int i = 0; i < A.get_rows(); i++) {
            for (int j = 0; j < A.get_columns(); j++) {
                Right->operator()(i, j) = (*this)(i, j);
            }
        }
    }

    void operator=(Matrix secMatrix) {
        for (int i = 0; i < secMatrix.get_rows(); i++) {
            for (int j = 0; j < secMatrix.get_columns(); j++) {
                (*this)(i, j) = secMatrix(i, j);
            }
        }
    }
};

class EliminationMatrix : public Matrix {
public:
    EliminationMatrix(int n, int i, int j, Matrix& A) : Matrix(n, n) {
        for (int k = 0; k < n; k++) {
            (*this)(k, k) = 1;
        }
        (*this)(i, j) = -A(i, j) / A(j, j);
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



Matrix inverse(Matrix matrix) {
    if (matrix.get_rows() != matrix.get_columns()) {
        cout << "Error: the dimensional problem occurred\n";
        return Matrix(0, 0);
    }
    if (matrix.determinant() == 0) {
        cout << "Error: the determinant is zero\n";
        return Matrix(0, 0);
    }
    Matrix result(matrix.get_size(), matrix.get_size());
    AugmentedMatrix augMatrix(matrix, IdentityMatrix(matrix.get_size()));
    for (int i = 0; i < matrix.get_size(); i++) {
        for (int j = 0; j < matrix.get_size(); j++) {
            if (i != j) {
                augMatrix = EliminationMatrix(matrix.get_size(), i, j, *augMatrix.Right) * augMatrix;
            }
        }
    }
    for (int i = 0; i < matrix.get_size(); i++) {
        for (int j = 0; j < matrix.get_size(); j++) {
            result(i, j) = augMatrix(i, j + matrix.get_size()) / augMatrix(i, i);
        }
    }
    return result;
}



#ifdef WIN32
    #define GNUPLOT_NAME "C:\\gnuplot\\bin\\gnuplot -persist"
#else
    #define GNUPLOT_NAME "gnuplot -persist"
#endif
int main() {
    #ifdef WIN32
    FILE* pipe = _popen(GNUPLOT_NAME, "w");
    #else
        FILE* pipe = popen(GNUPLOT_NAME, "w");
    #endif

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

    Matrix A_TA = A_transpose * A;

    Matrix A_inverse = inverse(A_TA);

    cout << "(A_T*A)^-1:\n" << A_inverse;

    Matrix A_transpose_B = A_transpose * B;

    cout << "A_T*b:\n" << A_transpose_B;

    Matrix x = A_inverse * A_transpose_B;

    cout << "x~:\n" << x;


    fprintf(pipe, "plot [-4 : 4] [-2 : 2] %lf*x**2 + %lf*x**1 + %lf*x**0 , '-' using 1:2 with points\n", x(0, 0), x(1, 0), x(2, 0));
    for (int i = 0; i < size; i++) {
        fprintf(pipe, "%d %d\n", t[i], b[i]);
    }
    fprintf(pipe, "e\n");
    fflush(pipe);
    

    #ifdef WIN32
        _pclose(pipe);
    #else
        pclose(pipe);
    #endif

    return 0;
}

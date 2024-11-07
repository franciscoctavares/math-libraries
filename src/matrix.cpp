#include "../include/matrix.h"

#include <iostream>
#include <iomanip>
#include <cmath>

Matrix::Matrix(Vector newElements, int a, int b) {
    if(a < 0 || b < 0) throw std::invalid_argument("Both arguments must be positive numbers");
    else if(a * b != newElements.size()) throw std::invalid_argument("Total number of elements(a * b) must correspond to the size of the vector");
    elements = newElements;
    n = a;
    m = b;
}

unsigned Matrix::getNRows(void) {
    return n;
}

unsigned Matrix::getNColumns(void) {
    return m;
}

void Matrix::printMatrix(void) {
    for(int i = 0; i < n; i++) {
        std::cout << "|";
        for(int j = 0; j < m; j++) {
            if(elements[i * m + j] < 0) std::cout << "-";
            else std::cout << "+";
            std::cout << std::setprecision(3) << std::fixed << fabs(elements[i * m + j]);
            if(j < m - 1) std::cout << " ";
        }
        std::cout << "|";
        if(i < n - 1) std::cout << std::endl;
    }
    std::cout << std::endl;
}

Matrix Matrix::operator+(Matrix matrix) {
    if(n != matrix.n || m != matrix.m) throw std::invalid_argument("Sizes of the 2 matrices must coincide");
    return Matrix(elements + matrix.elements, n, m);
}

void Matrix::operator+=(Matrix matrix) {
    if(n != matrix.n || m != matrix.m) throw std::invalid_argument("Sizes of the 2 matrices must coincide");
    elements += matrix.elements;
}

Matrix Matrix::operator*(double number) {
    return Matrix(elements * number, n, m);
}

void Matrix::operator*=(double number) {
    elements *= number;
}

Matrix Matrix::operator-(Matrix matrix) {
    if(n != matrix.n || m != matrix.m) throw std::invalid_argument("Sizes of the 2 matrices must coincide");
    return Matrix(elements - matrix.elements, n, m);
}

void Matrix::operator-=(Matrix matrix) {
    if(n != matrix.n || m != matrix.m) throw std::invalid_argument("Sizes of the 2 matrices must coincide");
    elements -= matrix.elements;
}

Vector Matrix::row(int row) {
    if(row < 0 || row >= n) throw std::invalid_argument("Argument must be positive and less than the number of rows");
    std::vector<double> aux;
    for(int i = 0; i < n; i++) {
        for(int j = 0; j < m; j++) {
            if(i == row) aux.push_back(elements[i * m + j]);
        }
    }
    return Vector(aux);
}

Vector Matrix::column(int column) {
    if(column < 0 || column >= m) throw std::invalid_argument("Argument must be positive and less than the number of columns");
    std::vector<double> aux;
    for(int i = 0; i < n; i++) {
        aux.push_back(elements[i * m + column]);
    }
    return Vector(aux);
}

Matrix Matrix::operator*(Matrix matrix) {
    if(m != matrix.n) throw std::invalid_argument("matrix sizes are not compatible");
    std::vector<double> newMatrix;
    for(int i = 0; i < n; i++) {
        for(int j = 0; j < matrix.m; j++) {
            newMatrix.push_back(row(i) * matrix.column(j));
        }
    }
    return Matrix(Vector(newMatrix), n, matrix.m);
}

void Matrix::operator*=(Matrix matrix) {
    elements = operator*(matrix).elements;
}

Matrix Matrix::transpose(void) {
    Vector newElements;
    for(int i = 0; i < m; i++) {
        newElements = newElements.attach(column(i));
    }
    return Matrix(newElements, m, n);
}

Matrix Matrix::stackVertical(Matrix matrix) {
    if(m != matrix.m) throw std::invalid_argument("Both matrices must have the same number of columns");
    return Matrix(elements.attach(matrix.elements), n + matrix.n, m);
}

Matrix Matrix::stackHorizontal(Matrix matrix) {
    if(n != matrix.n) throw std::invalid_argument("Both matrices must have the same number of rows");
    Vector newElements;
    for(int i = 0; i < n; i++) {
        newElements = newElements.attach(row(i).attach(matrix.row(i)));
    }
    return Matrix(newElements, n, m + matrix.m);
}

Matrix Matrix::subMatrix(unsigned a, unsigned b, unsigned c, unsigned d) {
    std::vector<double> newElements;
    if(!(a >= 0 && b < n && a <= b && c >= 0 && d < n && c <= d)) throw std::invalid_argument("Invalid indices");
    for(int i = 0; i < n; i++) {
        for(int j = 0; j < m; j++) {
            if(i >= a && i <= b && j >= c && j <= d) {
                newElements.push_back(elements[i * m + j]);
            }
        }
    }
    return Matrix(Vector(newElements), b - a + 1, d - c + 1);
}

Matrix Matrix::removeRow(int a) {
    // handle possible errors
    if(a < 0 || a >= n) throw std::invalid_argument("row must be a positive number and less than the number of rows");

    if(a == 0) return subMatrix(1, n - 1, 0, m - 1);
    else if(a == n - 1) return subMatrix(0, n - 2, 0, m - 1);
    else return subMatrix(0, a - 1, 0, m - 1).stackVertical(subMatrix(a + 1, n - 1, 0, m - 1));
}   

Matrix Matrix::removeColumn(int b) {
    //handle possible errors
    if(b < 0 || b >= m) throw std::invalid_argument("column must be a positive number and less than the number of columns");
    
    if(b == 0) return subMatrix(0, n - 1, 1, m - 1);
    else if(b == m - 1) return subMatrix(0, n - 1, 0, m - 2);
    else return subMatrix(0, n - 1, 0, b - 1).stackHorizontal(subMatrix(0, n - 1, b + 1, m - 1));
}

bool Matrix::isSquare(void) {
    return (n == m) ? true : false;
}

int Matrix::determinant(void) {
    if(!isSquare()) throw std::invalid_argument("The matrix must be a square matrix");
    else if(n < 2) throw std::invalid_argument("The matrix must be at least 2x2");

    if(n == 2) {
        return elements[0] * elements[3] - elements[1] * elements[2];
    }
    else {
        unsigned determinant = 0;

        for(int i = 0; i < m; i++) {
            determinant += pow(-1, i) * elements[0 * m + i] * removeColumn(i).removeRow(0).determinant();
        }

        return determinant;
    }
}

bool Matrix::isMatrixTall(void) {
    return (n > m) ? true : false;
}

bool Matrix::isMatrixWide(void) {
    return (m > n) ? true : false;
}

Matrix Matrix::rowOperation(int targetRow, int originRow, double factor) {
    std::vector<double> newElements;

    for(int i = 0; i < n; i++) {
        for(int j = 0; j < m; j++) {
            if(i != targetRow) {
                newElements.push_back(elements[i * m + j]);
            }
            else newElements.push_back(elements[targetRow * m + j] + factor * elements[originRow * m + j]);
        }
    }

    return Matrix(Vector(newElements), n, m);
}

Matrix Matrix::columnOperation(int targetColumn, int originColumn, double factor) {
    std::vector<double> newElements;

    for(int i = 0; i < n; i++) {
        for(int j = 0; j < m; j++) {
            if(j != targetColumn) newElements.push_back(elements[i * m + j]);
            else newElements.push_back(elements[i * m + targetColumn] + factor * elements[i * m + originColumn]);
        }
    }

    return Matrix(Vector(newElements), n, m);
}

Matrix Matrix::gaussElimination(void) {
    Matrix newMatrix(elements, n, m);
    //std::cout << elements[n] << " " << elements[0] << std::endl;
    //return newMatrix.rowOperation(1, 0, -1 * elements[m] * elements[0]);
    int originalRow = 0, targetRow = 0;

    //for(int i = 0; i < m; i++) {
    //    newMatrix = rowOperation();
    //    originalRow++;
    //}

    return newMatrix;
}

double Matrix::getElement(int row, int column) {
    if(row < 0 || row > n - 1 || column < 0 || column > m - 1) throw std::invalid_argument("row and column must both be positive numbers");
    return elements[row * m + column];
}

// static method
Matrix Matrix::identityMatrix(int n) {
    if(n <= 1) throw std::invalid_argument("Matrix size must be a positive number greater than 1");
    else {
        std::vector<double> newElements;
        for(int i = 0; i < n; i ++) {
            for(int j = 0; j < n; j++) {
                if(i == j) newElements.push_back(1);
                else newElements.push_back(0);
            }
        }
        return Matrix(Vector(newElements), n, n);
    }
}
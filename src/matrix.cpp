#include "../include/matrix.h"
#include <cmath>
#include <iostream>
#include <iomanip>
#include <algorithm>

Matrix zeros(unsigned rows, unsigned columns) {
    std::vector<double> vec;
    for(int i = 0; i < rows * columns; i++) {
        vec.push_back(0.0);
    }
    return Matrix(vec, rows, columns);
}

Matrix basisVector(unsigned size, unsigned index) {
    std::vector<double> vec;
    for(int i = 0; i < size; i++) {
        if(i == index) vec.push_back(1.0);
        else vec.push_back(0.0);
    }
    return Matrix(vec, size, 1);
}

Matrix::Matrix(std::vector<double> vec, unsigned a, unsigned b) {
    elements = vec;
    if(elements.size() != a * b) {
        n = elements.size();
        m = 1;
    }
    else {
        n = a;
        m = b;
    }
}

void Matrix::displayMatrix() {
    for(int i = 0; i < n; i++) {
        std::cout << "|";
        for(int j = 0; j < m; j++) {
            if(elements[i * m + j] < 0) std::cout << "-";
            else std::cout << "+";
            std::cout << std::setprecision(3) << std::fixed << fabs(elements[i * m + j]);
            if(j < m - 1) std::cout << " ";
        }
        std::cout << "|" << std::endl;
    }
}

Matrix Matrix::operator+(Matrix matrix) {
    std::ostringstream errorMsg;
    if(n != matrix.n || m != matrix.m) {
        errorMsg << "Error using operator+: matrix1 has dimensions (" << n << ", " << m << 
                    ") and matrix 2 has dimenions (" << matrix.n << ", " << matrix.m << ")";
        throw std::runtime_error(errorMsg.str());
    }
    Matrix aux(elements, n, m);
    for(int i = 0; i < n * m; i++) {
        aux.elements[i] += matrix.elements[i];
    }
    return aux;
}

Matrix Matrix::operator-(Matrix matrix) {
    std::ostringstream errorMsg;
    if(n != matrix.n || m != matrix.m) {
        errorMsg << "Error using operator-: matrix1 has dimensions (" << n << ", " << m << 
                    ") and matrix 2 has dimenions (" << matrix.n << ", " << matrix.m << ")";
        throw std::runtime_error(errorMsg.str());
    }
    Matrix aux(elements, n, m);
    for(int i = 0; i < n * m; i++) {
        aux.elements[i] -= matrix.elements[i];
    }
    return aux;
}

Matrix Matrix::operator*(Matrix matrix) {
    if(m == matrix.n) {
        std::vector<double> newStuff;
        for(int i = 0; i < n * matrix.m; i++) newStuff.push_back(0.0);
        Matrix newMatrix(newStuff, n, matrix.m);
        double aux;
        for(int i = 0; i < n; i++) {
            for(int j = 0; j < matrix.m; j++) {
                aux = 0;
                for(int k = 0; k < m; k++) {
                    aux += elements[i * m + k] * matrix.elements[k * matrix.n + j];
                }
                newMatrix.elements[i * matrix.m + j] = aux;
            }
        }
        return newMatrix;
    }
    else throw std::runtime_error("Dimensions don't match");
}

Matrix Matrix::getRow(unsigned r) {
    if(r >= 0 && r < n) {
        std::vector<double> aux;

        for(int j = 0; j < m; j++) {
            aux.push_back(elements[r * m + j]);
        }
        return Matrix(aux, 1, m);
    }
    else {
        std::ostringstream errorMsg;
        errorMsg << "Error using getRow: the argument must be between 0 and " << n - 1 << ", but the value provided was " << r;
        throw std::runtime_error(errorMsg.str());
    }
}

Matrix Matrix::getColumn(unsigned c) {
    if(c >= 0 && c < m) {
        std::vector<double> aux;

        for(int i = 0; i < n; i++) {
            aux.push_back(elements[i * m + c]);
        }
        return Matrix(aux, n, 1);
    }
    else {
        std::ostringstream errorMsg;
        errorMsg << "Error using getColumn: the argument must be between 0 and " << m - 1 << ", but the value provided was " << c;
        throw std::runtime_error(errorMsg.str());
    }
}

Matrix Matrix::rowOperation(unsigned a, unsigned b, double factor) {
    Matrix newMatrix(elements, n, m);
    // a = source, b = target
    for(int j = 0; j < m; j++) {
        newMatrix.elements[b * m + j] = newMatrix.elements[b * m + j] + newMatrix.elements[a * m + j] * factor;
    }

    return newMatrix;
}

Matrix Matrix::columnOperation(unsigned a, unsigned b, double factor) {
    Matrix newMatrix(elements, n, m);
    // a = source, b = target
    for(int i = 0; i < n; i++) {
        newMatrix.elements[i * m + b] = newMatrix.elements[i * m + b] + newMatrix.elements[i * m + a] * factor;
    }

    return newMatrix;
}

double Matrix::dotProduct(Matrix matrix) { // assumes both matrices are row/column matrices(matrix 1 could be row and 2 could be column)
    if((n != 1 && m != 1) || (matrix.n != 1 && matrix.m != 1)) {
        std::ostringstream errorMsg;
        errorMsg << "Error using dotProduct: Both matrices have to be vector(row or column) matrices";
        throw std::runtime_error(errorMsg.str());
    }
    else if(elements.size() != matrix.elements.size()) {
        std::ostringstream errorMsg;
        errorMsg << "Error using dotProduct: cannot calculate the dot product when matrix 1 and matrix 2 have " << elements.size() 
                 << " and " << matrix.elements.size() << " elements, respectively";
        throw std::runtime_error(errorMsg.str());
    }

    double result = 0;
    for(int i = 0; i < elements.size(); i++) {
        result += elements[i] * matrix.elements[i];
    }
    return result;

}

unsigned Matrix::rows() {
    return n;
}

unsigned Matrix::columns() {
    return m;
}

double Matrix::getElement(unsigned row, unsigned column) {
    return elements[row * m + column];
}

void Matrix::setElement(unsigned row, unsigned col, double value) {
    elements[row * m + col] = value;
}

Matrix Matrix::stackVertical(Matrix matrix) {
    if(m != matrix.m) {
        std::ostringstream errorMsg;
        errorMsg << "Error using stackVertical: cannot vertically stack a matrix with " << matrix.m << " columns below a matrix with " << m << " columns";
        throw std::runtime_error(errorMsg.str());
    }
    std::vector<double> aux;
    for(int i = 0; i < rows() * columns(); i++) aux.push_back(getElement(i / rows(), i % columns()));
    for(int j = 0; j < matrix.rows() * matrix.columns(); j++) aux.push_back(matrix.getElement(j / matrix.rows(), j % matrix.columns()));
    return Matrix(aux, n + matrix.rows(), m);
}

Matrix Matrix::stackHorizontal(Matrix matrix) {
    if(n != matrix.n) {
        std::ostringstream errorMsg;
        errorMsg << "Error using stackHorizontal: cannot horizontally stack a matrix with " << matrix.n << " rows next to a matrix with " << n << " rows";
        throw std::runtime_error(errorMsg.str());
    }
    std::vector<double> aux;
    for(int i = 0; i < n; i++) {
        for(int j = 0; j < m; j++) {
            aux.push_back(elements[i * m + j]);
        }
        for(int k = 0; k < matrix.columns(); k++) {
            aux.push_back(matrix.getElement(i, k));
        }
    }
    return Matrix(aux, n, m + matrix.columns());
}

Matrix Matrix::transpose() {
    Matrix newMatrix = zeros(m, n);

    for(int i = 0; i < n; i++) {
        for(int j = 0; j < m; j++) {
            newMatrix.setElement(j, i, getElement(i, j));
        }
    }
    return newMatrix;
}

Matrix Matrix::operator*(double value) {
    std::vector<double> aux = elements;
    for(int i = 0; i < n * m; i++) aux[i] *= value;
    return Matrix(aux, n, m);
}

unsigned Matrix::maxValueIndex() {
    auto maxIt = std::max_element(elements.begin(), elements.end());
    size_t maxIndex = std::distance(elements.begin(), maxIt);
    return maxIndex;
}

unsigned Matrix::minValueIndex() {
    auto minIt = std::min_element(elements.begin(), elements.end());
    size_t minIndex = std::distance(elements.begin(), minIt);
    return minIndex;
}

Matrix Matrix::pointDivision(Matrix matrix) {
    if(rows() != matrix.rows() || columns() != matrix.columns()) {
        std::ostringstream errorMsg;
        errorMsg << "Error using pointDivision: cannot multiply matrix1(" << n << " x " << m << ") by matrix2(" << matrix.n << " x " << matrix.m << ")";
        throw std::runtime_error(errorMsg.str());
    }
    std::vector<double> aux = elements;
    for(int i = 0; i < rows() * columns(); i++) {
        aux[i] /= matrix.getElement(i / m, i % m);
    }
    return Matrix(aux, n, m);
}

Matrix Matrix::setRow(unsigned row, Matrix matrix) {
    if(row < 0 || row > rows() - 1) {
        std::ostringstream errorMsg;
        errorMsg << "Error using setRow: the row argument must be between 0 and " << n - 1 << ", but the value provided was " << row;
        throw std::runtime_error(errorMsg.str());
    }
    else if(columns() != matrix.columns()) {
        std::ostringstream errorMsg;
        errorMsg << "Error using setRow: the matrix has " << columns() << " columns, but the new row matrix has " << matrix.columns() << " columns";
        throw std::runtime_error(errorMsg.str());
    }

    std::vector<double> aux;
    for(int i = 0; i < n; i++) {
        for(int j = 0; j < m; j++) {
            if(i == row) aux.push_back(matrix.getElement(0, j));
            else aux.push_back(elements[i * m + j]);
        }
    }
    return Matrix(aux, n, m);
}

Matrix Matrix::setColumn(unsigned column, Matrix matrix) {
    if(column < 0 || column > columns() - 1) {
        std::ostringstream errorMsg;
        errorMsg << "Error using setColumn: the column argument must be between 0 and " << m - 1 << ", but the value provided was " << column;
        throw std::runtime_error(errorMsg.str());
    }
    else if(rows() != matrix.rows()) {
        std::ostringstream errorMsg;
        errorMsg << "Error using setRow: the matrix has " << rows() << " rows, but the new column matrix has " << matrix.rows() << " rows";
        throw std::runtime_error(errorMsg.str());
    }
    
    std::vector<double> aux;
    for(int i = 0; i < n; i++) {
        for(int j = 0; j < m; j++) {
            if(j == column) aux.push_back(matrix.getElement(i, 0));
            else aux.push_back(elements[i * m + j]);
        }
    }
    return Matrix(aux, n, m);
}
#include "../include/matrix.h"
#include <cmath>
#include <iostream>
#include <iomanip>
#include <algorithm>

#include <string>

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
    size_t maxWidth = 0;
    for(double val: elements) {
        std::string str = std::to_string(val);
        // Trim trailing zeroes for nicer formatting (optional)
        str.erase(str.find_last_not_of('0') + 1, std::string::npos);
        if (str.back() == '.') str.pop_back(); // remove trailing dot if needed
        maxWidth = std::max(maxWidth, str.length());
    }

    //maxWidth += 2;

    for(int i = 0; i < n; i++) {
        std::cout << "|";
        for(int j = 0; j < m; j++) {
            //if(elements[i * m + j] < 0) std::cout << "-";
            //else std::cout << "+";
            std::cout << std::setw(maxWidth) << elements[i * m + j];
            //std::cout << std::setprecision(3) << std::fixed << fabs(elements[i * m + j]);
            //if(j < m - 1) std::cout << " ";
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

void Matrix::operator+=(Matrix matrix) {
    std::ostringstream errorMsg;
    if(n != matrix.n || m != matrix.m) {
        errorMsg << "Error using operator+=: matrix1 has dimensions (" << n << ", " << m << 
                    ") and matrix 2 has dimenions (" << matrix.n << ", " << matrix.m << ")";
        throw std::runtime_error(errorMsg.str());
    }
    for(int i = 0; i < n; i++) {
        for(int j = 0; j < m; j++) {
            elements[i * m + j] += matrix.elements[i * m + j];
        }
    }
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

void Matrix::operator-=(Matrix matrix) {
    std::ostringstream errorMsg;
    if(n != matrix.n || m != matrix.m) {
        errorMsg << "Error using operator-=: matrix1 has dimensions (" << n << ", " << m << 
                    ") and matrix 2 has dimenions (" << matrix.n << ", " << matrix.m << ")";
        throw std::runtime_error(errorMsg.str());
    }
    for(int i = 0; i < n; i++) {
        for(int j = 0; j < m; j++) {
            elements[i * m + j] -= matrix.elements[i * m + j];
        }
    }
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
/*
void Matrix::operator=(Matrix matrix) {
    elements.clear();
    n = matrix.n;
    m = matrix.m;
    for(int i = 0; i < n * m; i++) {
        elements.push_back(matrix.elements[i]);
    }
}
*/


std::vector<double> Matrix::getRow(unsigned r) {
    if(r >= 0 && r < n) {
        std::vector<double> aux;

        for(int j = 0; j < m; j++) {
            aux.push_back(elements[r * m + j]);
        }
        return aux;
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

void Matrix::rowOperation(unsigned a, unsigned b, double factor) {
    if(a < 0 || a >= n) {
        std::ostringstream errorMsg;
        errorMsg << "Error using rowOperation: The a argument must be between 0 and " << n - 1 << ", but the value provided was " << a;
        throw std::runtime_error(errorMsg.str());
    }
    if(b < 0 || b >= n) {
        std::ostringstream errorMsg;
        errorMsg << "Error using rowOperation: The b argument must be between 0 and " << n - 1 << ", but the value provided was " << b;
        throw std::runtime_error(errorMsg.str());
    }
    for(int j = 0; j < m; j++) {
        elements[b * m + j] += elements[a * m + j] * factor;
    }
}

void Matrix::columnOperation(unsigned a, unsigned b, double factor) {
    /*
    Matrix newMatrix(elements, n, m);
    // a = source, b = target
    for(int i = 0; i < n; i++) {
        newMatrix.elements[i * m + b] = newMatrix.elements[i * m + b] + newMatrix.elements[i * m + a] * factor;
    }

    return newMatrix;
    */

    for(int i = 0; i < n; i++) {
        elements[i * m + b] += elements[i * m + a] * factor;
    }
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
    if(row < 0 || row >= n) {
        std::ostringstream errorMsg;
        errorMsg << "Error using setElement: The row argument must be between 0 and " << n - 1 << ", but the value provided was " << row;
        throw std::runtime_error(errorMsg.str());
    }
    if(col < 0 || col >= m) {
        std::ostringstream errorMsg;
        errorMsg << "Error using setElement: The col argument must be between 0 and " << m - 1 << ", but the value provided was " << col;
        throw std::runtime_error(errorMsg.str());
    }
    elements[row * m + col] = value;
}

void Matrix::stackVertical(Matrix matrix) {
    if(m != matrix.m) {
        std::ostringstream errorMsg;
        errorMsg << "Error using stackVertical: cannot vertically stack a matrix with " << matrix.m << " columns below a matrix with " << m << " columns";
        throw std::runtime_error(errorMsg.str());
    }
    /*
    std::vector<double> aux;
    for(int i = 0; i < rows() * columns(); i++) aux.push_back(getElement(i / rows(), i % columns()));
    for(int j = 0; j < matrix.rows() * matrix.columns(); j++) aux.push_back(matrix.getElement(j / matrix.rows(), j % matrix.columns()));
    return Matrix(aux, n + matrix.rows(), m);
    */
    for(int i = 0; i < matrix.n * matrix.m; i++) {
        elements.push_back(matrix.elements[i]);
    }
    n += matrix.n;
}

void Matrix::stackHorizontal(Matrix matrix) {
    if(n != matrix.n) {
        std::ostringstream errorMsg;
        errorMsg << "Error using stackHorizontal: cannot horizontally stack a matrix with " << matrix.n << " rows next to a matrix with " << n << " rows";
        throw std::runtime_error(errorMsg.str());
    }

    unsigned current_index;
    for(int i = 0; i < n; i++) {
        //std::cout << elements.size() << std::endl;
        for(int j = 0; j < matrix.m; j++) {
            current_index = i * (m + matrix.m) + (m + j);
            elements.insert(elements.begin() + current_index, matrix.elements[i * matrix.m + j]);
            //std::cout << current_index << " ";
        }
        //std::cout << std::endl;
    }
    m += matrix.m;
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

void Matrix::operator*=(double value) {
    for(int i = 0; i < n * m; i++) elements[i] *= value;
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

Matrix Matrix::removeRow(unsigned row) {
    if(row < 0 || row > rows() - 1) {
        std::ostringstream errorMsg;
        errorMsg << "Error using removeRow: the row argument must be between 0 and " << n - 1 << ", but the value provided was " << row;
        throw std::runtime_error(errorMsg.str());
    }

    std::vector<double> aux;
    for(int i = 0; i < rows(); i++) {
        if(i == row) continue;
        for(int j = 0; j < columns(); j++) {
            aux.push_back(getElement(i, j));
        }
    }
    return Matrix(aux, n - 1, m);
}

void Matrix::removeColumn(unsigned column) {
    if(column < 0 || column > columns() - 1) {
        std::ostringstream errorMsg;
        errorMsg << "Error using removeColumn: the column argument must be between 0 and " << m - 1 << ", but the value provided was " << column;
        throw std::runtime_error(errorMsg.str());
    }

    std::vector<double> aux = elements;
    elements.clear();
    for(int i = 0; i < rows(); i++) {
        for(int j = 0; j < columns(); j++) {
            if(j == column) continue;
            elements.push_back(getElement(i, j));
        }
    }
    //return Matrix(aux, n, m - 1);
    m--;
}

unsigned Matrix::findValueInVectorMatrix(double value) {
    for(int i = 0; i < elements.size(); i++) {
        if(elements[i] == value) return i;
    }
    return -1;
}

Matrix Matrix::subMatrix(unsigned r1, unsigned r2, unsigned c1, unsigned c2) {
    std::vector<double> aux;

    for(int i = 0; i < n; i++) {
        for(int j = 0; j < m; j++) {
            if(i >= r1 && i <= r2 && j >= c1 && j <= c2) aux.push_back(elements[i * m + j]);
        }
    }
    return Matrix(aux, r2 - r1 + 1, c2 - c1 + 1);
}

std::vector<double> Matrix::getElements() {
    return elements;
}
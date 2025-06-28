#include "../include/matrix.h"
#include <cmath>
#include <iostream>
#include <iomanip>

Matrix zeros(unsigned rows, unsigned columns) {
    std::vector<double> vec;
    for(int i = 0; i < rows * columns; i++) {
        vec.push_back(0.0);
    }
    return Matrix(vec, rows, columns);
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
    std::vector<double> v = {0};
    Matrix aux(elements, n, m);
    if(n != matrix.n || m != matrix.m) return Matrix(v, 1, 1);
    for(int i = 0; i < n * m; i++) {
        aux.elements[i] += matrix.elements[i];
    }
    return aux;
}

Matrix Matrix::operator-(Matrix matrix) {
    std::vector<double> v = {0};
    Matrix aux(elements, n, m);
    if(n != matrix.n || m != matrix.m) return Matrix(v, 1, 1);
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
}

Matrix Matrix::getRow(unsigned r) {
    if(r >= 0 && r < n) {
        std::vector<double> aux;

        for(int j = 0; j < m; j++) {
            aux.push_back(elements[r * m + j]);
        }
        return Matrix(aux, 1, m);
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

double Matrix::dotProduct(Matrix matrix) {
    if(elements.size() == matrix.elements.size()) {
        double result = 0;
        for(int i = 0; i < elements.size(); i++) {
            result += elements[i] * matrix.elements[i];
        }
        return result;
    }
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

Matrix Matrix::stackVertical(Matrix matrix) {
    if(m == matrix.columns()) {
        std::vector<double> aux;
        for(int i = 0; i < rows() * columns(); i++) aux.push_back(getElement(i / rows(), i % columns()));
        for(int j = 0; j < matrix.rows() * matrix.columns(); j++) aux.push_back(matrix.getElement(j / matrix.rows(), j % matrix.columns()));
        return Matrix(aux, n + matrix.rows(), m);
    }
}

Matrix Matrix::stackHorizontal(Matrix matrix) {
    if(n == matrix.rows()) {
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
}
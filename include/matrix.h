#ifndef MATRIX_H
#define MATRIX_H

#include <vector>

class Matrix {
    private:
        std::vector<double> elements;
        unsigned n, m;
    public:
        Matrix(std::vector<double>, unsigned, unsigned);
        void displayMatrix();
        Matrix operator+(Matrix);
        Matrix operator-(Matrix);
        Matrix operator*(Matrix);
        Matrix getRow(unsigned);
        Matrix getColumn(unsigned);
        Matrix rowOperation(unsigned, unsigned, double);
        Matrix columnOperation(unsigned, unsigned, double);
};

#endif
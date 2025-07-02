#ifndef MATRIX_H
#define MATRIX_H

#include <vector>

class Matrix {
    private:
        std::vector<double> elements;
        unsigned n, m;
    public:
        Matrix(void) = default;
        Matrix(std::vector<double>, unsigned, unsigned);
        void displayMatrix();
        Matrix operator+(Matrix);
        Matrix operator-(Matrix);
        Matrix operator*(Matrix);
        Matrix getRow(unsigned);
        Matrix getColumn(unsigned);
        Matrix rowOperation(unsigned, unsigned, double);
        Matrix columnOperation(unsigned, unsigned, double);
        double dotProduct(Matrix); // dot product, when the 2 matrices are vectors
        unsigned rows(); // returns the number of rows
        unsigned columns(); // returns the number of columns
        double getElement(unsigned, unsigned); // returns the element given the indices of the row and column
        void setElement(unsigned, unsigned, double); // given the indices of the row and column, sets the element to the value specified by the argument
        Matrix stackVertical(Matrix);
        Matrix stackHorizontal(Matrix);
        Matrix transpose();
        Matrix operator*(double);
        unsigned maxValueIndex();
        unsigned minValueIndex();
        Matrix pointDivision(Matrix);
        Matrix setRow(unsigned, Matrix);
        Matrix setColumn(unsigned, Matrix);
};

Matrix zeros(unsigned, unsigned);
Matrix basisVector(unsigned, unsigned);

#endif
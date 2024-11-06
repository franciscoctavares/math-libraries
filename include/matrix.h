#ifndef MATRIX_H
#define MATRIX_H

#include "vector.h"

class Matrix {
    private:
        Vector elements;
        unsigned n, m;
    public:
        Matrix(Vector, int, int);
        /**
         * Returns the number of rows
        */
        unsigned getNRows(void);
        /**
         * Returns the number of columns
        */
        unsigned getNColumns(void);
        /**
         * Prints the matrix to the terminal
        */
        void printMatrix(void);
        /**
         * Adds 2 matrices, returning the result
        */
        Matrix operator+(Matrix);
        /**
         * Adds 2 matrices, and the result is stored in the caller object
        */
        void operator+=(Matrix);
        // multiply matrix by number
        Matrix operator*(double);
        void operator*=(double);
        // subtract 2 matrices
        Matrix operator-(Matrix);
        void operator-=(Matrix);
        // obtain a row or a column from a matrix
        Vector row(int);
        Vector column(int);
        // multiply 2 matrices
        Matrix operator*(Matrix);
        void operator*=(Matrix);
        /**
         * Returns the transpose of the matrix
        */
        Matrix transpose(void);
        /**
         * Stacks 2 matrices vertically
        */
        Matrix stackVertical(Matrix);
        /**
         * Stacks 2 matrices horizontally
        */
        Matrix stackHorizontal(Matrix);
        Matrix subMatrix(unsigned, unsigned, unsigned, unsigned);
        /**
         * Removes the specified row from the matrix
        */
        Matrix removeRow(int);
        /**
         * Removes the specified column from the matrix
        */
        Matrix removeColumn(int);
        /**
         * Simple helper method that returns true if the matrix is a square matrix and false otherwise.
        */
        bool isSquare(void);
        /** 
         *  Calculates the determinant of the matrix.
         * @exception std::invalid_argument, if the matrix is not square
         */
        int determinant(void);
        bool isMatrixTall(void);
        bool isMatrixWide(void);
        Matrix rowOperation(int, int, double);
        Matrix columnOperation(int, int, double);
        Matrix gaussElimination(void);
        double getElement(int, int);
        
        static Matrix identityMatrix(int);
};

#endif
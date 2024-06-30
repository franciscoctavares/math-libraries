#ifndef MATRIX_H
#define MATRIX_H

#include "../vector/vector.h"

class Matrix {
    private:
        Vector elements;
        unsigned int n, m;
    public:
        Matrix(Vector, unsigned int, unsigned int);
        void printMatrix(void);
};

#endif
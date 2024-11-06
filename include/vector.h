#ifndef VECTOR_H
#define VECTOR_H

#include <vector>

class Vector {
    private:
        std::vector<double> elements;
    public:
        Vector() = default; // no need to test this
        Vector(std::vector<double>);
        Vector(int);
        unsigned size(void);
        void printVector(void); // no need to test this
        double operator[](int);
        // BEGIN OPERATORS
        Vector operator+(Vector);
        void operator+=(Vector);
        Vector operator*(double);
        void operator*=(double);
        Vector operator-(Vector);
        void operator-=(Vector);
        // END OPERATORS
        Vector attach(Vector);
        double operator*(Vector);
        Vector subVector(int, int);
};

#endif
#ifndef VECTOR_H
#define VECTOR_H

#include <vector>

using namespace std;

class Vector {
    private:
        vector<double> elements;
    public:
        Vector() = default;
        Vector(vector<double>);
        void printVector(void);
        Vector operator+(Vector);
        Vector operator-(Vector);
        Vector operator*(double);
        double operator*(Vector);
        Vector subVector(unsigned int, unsigned int);
        unsigned int size(void);
        double operator[](unsigned int);
};

#endif
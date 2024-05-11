#ifndef VECTOR_H
#define VECTOR_H

#include <vector>
#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;

class Vector {
    private:
        vector<double> elements;
    public:
        Vector(vector<double>);
        void printVector(void);
};

#endif
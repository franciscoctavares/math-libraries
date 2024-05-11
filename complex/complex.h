#ifndef COMPLEX_H
#define COMPLEX_H

#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;

class Complex {
    private:
        double real, imaginary;
    public:
        Complex(double a = 0.0, double b = 0.0);
        void printNumber(void);
        Complex operator+(Complex);
        Complex operator-(Complex);
        Complex operator*(Complex);
        Complex operator*(double);
        Complex operator/(Complex);
        Complex operator/(double);
};

#endif
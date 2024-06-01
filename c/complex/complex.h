#ifndef COMPLEX_H
#define COMPLEX_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

typedef struct {
    double real, imaginary;
}Complex;

void createNumber(Complex * number, double a, double b);
void printNumber(Complex number);
void addNumbers(Complex * result, Complex number1, Complex number2);
void subtractNumbers(Complex * result, Complex number1, Complex number2);

#endif
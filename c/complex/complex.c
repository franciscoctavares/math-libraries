#include "complex.h"

void createNumber(Complex * number, double a, double b) {
    number -> real= a;
    number -> imaginary = b;
}

void printNumber(Complex number) {
    char sign_a = (number.real < 0) ? '-' : '+';
    char sign_b = (number.imaginary < 0) ? '-' : '+';
    printf("%c%.3f %c%.3fi", sign_a, fabs(number.real), sign_b, fabs(number.imaginary));
}

void addNumbers(Complex * result, Complex number1, Complex number2) {
    result -> real = number1.real + number2.real;
    result -> imaginary = number1.imaginary + number2.imaginary;
}

void subtractNumbers(Complex * result, Complex number1, Complex number2) {
    result -> real = number1.real - number2.real;
    result -> imaginary = number1.imaginary - number2.imaginary;
}
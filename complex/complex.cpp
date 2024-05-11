#include "complex.h"

Complex :: Complex(double a, double b) {
    real = a;
    imaginary = b;
}

void Complex :: printNumber(void) {
    if(real < 0) cout << "-";
    else cout << "+";
    cout << " " << fixed << setprecision(3) << fabs(real) << " ";
    if(imaginary < 0) cout << "-";
    else cout << "+";
    cout << " " << fabs(imaginary) << "i";
}
    
Complex Complex :: operator+(Complex number) {
    double a = real, b = imaginary;
    double c = number.real, d = number.imaginary;
    return Complex(a + c, b + d);
}

Complex Complex :: operator-(Complex number) {
    double a = real, b = imaginary;
    double c = number.real, d = number.imaginary;
    return Complex(a - c, b - d);
}

Complex Complex :: operator*(Complex number) {
    double a = real, b = imaginary;
    double c = number.real, d = number.imaginary;
    return Complex(a * c - b * d, a * d + b* c);
}

Complex Complex ::  operator*(double number) {
    double a = real, b = imaginary;
    return Complex(a * number, b * number);
}

Complex Complex :: operator/(Complex number) {
    double a = real, b = imaginary;
    double c = number.real, d = number.imaginary;
    double denominator = c * c + d * d;
    return Complex((a * c + b * d) / denominator, (c * b - d * a) / denominator);
}

Complex Complex :: operator/(double number) {
    double a = real, b = imaginary;
    return Complex(a / number, b / number);
}
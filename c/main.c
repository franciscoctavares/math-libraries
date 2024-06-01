#include <stdio.h>
#include "complex/complex.h"

int main(int argc, char ** argv) {
    Complex num, num2;
    createNumber(&num, -1.4, -2.8);
    createNumber(&num2, 1.0, 2.0);
    subtractNumbers(&num2, num, num2);
    printNumber(num2);
    return 0;
}
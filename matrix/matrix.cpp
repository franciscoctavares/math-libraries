#include "matrix.h"
#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;

Matrix :: Matrix(Vector vector, unsigned int a, unsigned int b) {
    elements = vector;
    if(a * b != elements.size()) {
        n = elements.size();
        m = 1;
    }
    else {
        n = a;
        m = b;
    }
}

void Matrix :: printMatrix(void) {
    for(int i = 0; i < n; i++) {
        cout << "|";
        for(int j = 0; j < m; j++) {
            if(elements[i * m + j] < 0) cout << "-";
            else cout << "+";
            cout << setprecision(3) << fixed << fabs(elements[i * m + j]);
            if(j < m - 1) cout << " ";
        }
        cout << "|";
        if(i < n - 1) cout << endl;
    }
}
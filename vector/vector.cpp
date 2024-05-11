#include "vector.h"

Vector :: Vector(vector<double> newElements) {
    elements = newElements;
}

void Vector :: printVector(void) {
    cout << "[";
    for(int i = 0; i < elements.size(); i++) {
        if(elements[i] < 0) cout << "-";
        else cout << "+";
        cout << fixed << setprecision(3) << fabs(elements[i]);
        if(i < elements.size() - 1) cout << " ";
    }
    cout << "]";
}
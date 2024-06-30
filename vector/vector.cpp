#include "vector.h"
#include <iostream>
#include <iomanip>
#include <cmath>

Vector :: Vector(vector<double> newElements) {
    elements = newElements;
}

void Vector :: printVector(void) {
    cout << "[";
    for(int i = 0; i < elements.size(); i++) {
        if(elements[i] < 0) cout << "-";
        else cout << "+";
        cout << setprecision(3) << fixed << fabs(elements[i]);
        if(i < elements.size() - 1) cout << " ";
    }
    cout << "]";
}

Vector Vector :: operator+(Vector thing) {
    vector<double> newElements;
    if(elements.size() < thing.elements.size()) {
        for(int i = 0; i < elements.size(); i++) {
            newElements.push_back(elements[i] + thing.elements[i]);
        }
        for(int j = elements.size(); j < thing.elements.size(); j++) {
            newElements.push_back(thing.elements[j]);
        }
    }
    else {
        for(int i = 0; i < thing.elements.size(); i++) {
            newElements.push_back(thing.elements[i] + elements[i]);
        }
        for(int j = thing.elements.size(); j < elements.size(); j++) {
            newElements.push_back(elements[j]);
        }
    }
    return Vector(newElements);
}

Vector Vector :: operator-(Vector vector) {
    return operator+(vector * -1);
}

Vector Vector :: operator*(double number) {
    vector<double> newElements;
    for(int i = 0; i < elements.size(); i++) {
        newElements.push_back(elements[i] * number);
    }
    return newElements;
}

double Vector :: operator*(Vector vector) {
    double result = 0;
    if(elements.size() < vector.elements.size()) {
        for(int i = 0; i < elements.size(); i++) {
            result += elements[i] * vector.elements[i];
        }
    }
    else {
        for(int i = 0; i < vector.elements.size(); i++) {
            result += elements[i] * vector.elements[i];
        }
    }
    return result;
}

Vector Vector :: subVector(unsigned int a, unsigned int b) {
    vector<double> newElements;
    if(a < b && a >= 0 && a < elements.size() && b >= 0 && b < elements.size()) {
        for(int i = 0; i < elements.size(); i++) {
            if(i >= a && i <= b) newElements.push_back(elements[i]);
        }
    }
    else return Vector(elements);
    return Vector(newElements);
}

unsigned int Vector :: size(void) {
    return elements.size();
}

double Vector :: operator[](unsigned int index) {
    if(index >= 0 && index < elements.size()) return elements[index];
    else return NULL;
}
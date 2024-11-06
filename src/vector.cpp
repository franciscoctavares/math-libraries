#include "../include/vector.h"

#include <iostream>
#include <iomanip>
#include <cmath>

Vector :: Vector(std::vector<double> newElements) {
    elements = newElements;
}

Vector :: Vector(int nElements) {
    if(nElements < 0) throw std::invalid_argument("Size of vector must be a positive number");
    for(int i = 0; i < nElements; i++) elements.push_back(0.0);
}

unsigned Vector :: size(void) {
    return elements.size();
}

void Vector :: printVector(void) {
    std::cout << "[";
    for(int i = 0; i < elements.size(); i++) {
        if(elements[i] < 0) std::cout << "-";
        else std::cout << "+";
        std::cout << std::setprecision(3) << std::fixed << fabs(elements[i]);
        if(i < elements.size() - 1) std::cout << " ";
    }
    std::cout << "]" << std::endl;
}

double Vector :: operator[](int index) {
    //if(index < 0 || index >= elements.size()) throw std::invalid_argument("Index out of bounds");
    if(index < 0) throw std::invalid_argument("index must be a positive number");
    else if(index >= elements.size()) throw std::invalid_argument("index must be less than the vector's size");
    return elements[index];
}

Vector Vector :: operator+(Vector vectorOperator) {
    std::vector<double> newElements;
    if(elements.size() < vectorOperator.elements.size()) {
        for(int i = 0; i < elements.size(); i++) {
            newElements.push_back(elements[i] + vectorOperator.elements[i]);
        }
        for(int j = elements.size(); j < vectorOperator.elements.size(); j++) {
            newElements.push_back(vectorOperator.elements[j]);
        }
    }
    else {
        for(int i = 0; i < vectorOperator.elements.size(); i++) {
            newElements.push_back(elements[i] + vectorOperator.elements[i]);
        }
        for(int j = vectorOperator.elements.size(); j < elements.size(); j++) {
            newElements.push_back(elements[j]);
        }
    }
    return Vector(newElements);
}

void Vector :: operator+=(Vector vectorOperator) {
    elements = operator+(vectorOperator).elements;
}

Vector Vector :: operator*(double number) {
    std::vector<double> newElements;
    for(int i = 0; i < elements.size(); i++) {
        newElements.push_back(elements[i] * number);
    }
    return Vector(newElements);
}

void Vector :: operator*=(double number) {
    elements = Vector(operator*(number)).elements;
}

Vector Vector :: operator-(Vector vectorOperator) {
    return operator+(vectorOperator * -1);
}

void Vector :: operator-=(Vector vectorOperator) {
    elements = Vector(operator-(vectorOperator)).elements;
}

Vector Vector :: attach(Vector vectorOperator) {
    std::vector<double> newElements = elements;
    for(int i = 0; i < vectorOperator.elements.size(); i++) {
        newElements.push_back(vectorOperator.elements[i]);
    }
    return Vector(newElements);
}

double Vector :: operator*(Vector vectorOperator) {
    double sum = 0;
    if(size() < vectorOperator.size()) {
        for(int i = 0; i < size(); i++) {
            sum += (elements[i] * vectorOperator.elements[i]);
        }
    }
    else {
        for(int i = 0; i < vectorOperator.size(); i++) {
            sum += (elements[i] * vectorOperator.elements[i]);
        }
    }
    return sum;
}

Vector Vector :: subVector(int a, int b) {
    if(a < 0 || b < 0 || a > b || a >= elements.size() || b >= elements.size()) throw std::invalid_argument("One or both indices are out of bounds");
    std::vector<double> newElements;
    for(int i = a; i <= b; i++) {
        newElements.push_back(elements[i]);
    }
    return Vector(newElements);
}
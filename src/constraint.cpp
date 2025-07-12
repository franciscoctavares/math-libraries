#include "../include/constraint.h"

Constraint::Constraint(std::vector<double> newLhs, std::string newConstraintType, double newRhs) {
    lhs = newLhs;
    if(newConstraintType == "<=") type = LESS_THAN_OR_EQUAL;
    else if(newConstraintType == "=") type = EQUAL;
    else if(newConstraintType == ">=") type = GREATER_THAN_OR_EQUAL;
    else throw std::runtime_error("Invalid type of constraint. A constraint can only be of the 3 following types: <=, >= or =");
    rhs = newRhs;
}

std::vector<double> Constraint::getLhs() {
    return lhs;
}

ConstraintType Constraint::getType() {
    return type;
}

double Constraint::getRhs() {
    return rhs;
}
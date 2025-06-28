#ifndef LP_H
#define LP_H

#include "matrix.h"
#include <vector>

enum restrictionType {
    LESS_THAN_OR_EQUAL,      // <=
    EQUAL,                   // =
    GREATER_THAN_OR_EQUAL    // >=
};

class LpProblem {
    private:
        Matrix objectiveFunction; // row Matrix with the objective function coefficients
        Matrix restrictionsLHS, restrictionsRHS;
        std::vector<restrictionType> restrictionsTypes;
        Matrix solution;
    public:
        LpProblem(Matrix, Matrix, Matrix, std::vector<restrictionType>);
        bool isRestrictionSatisfied(Matrix, Matrix, double, restrictionType);
        bool isSolutionAdmissible(Matrix);
};

#endif
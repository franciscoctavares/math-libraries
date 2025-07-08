#ifndef LP_H
#define LP_H

#include "matrix.h"
#include <vector>

#include <limits>

#define M 1000000

enum restrictionType {
    LESS_THAN_OR_EQUAL,      // <=
    EQUAL,                   // =
    GREATER_THAN_OR_EQUAL    // >=
};

enum extraVariables {
    SLACK,                  // +s
    SURPLUS,                // -s
    NEUTRAL                 // +0
};

enum ProblemType {
    MIN,
    MAX
};

class LpProblem {
    private:
        Matrix objectiveFunction; // row Matrix with the objective function coefficients
        Matrix restrictionsLHS, restrictionsRHS;
        std::vector<restrictionType> restrictionsTypes;
        Matrix optimalSolution;
        ProblemType type;
        // private methods
        bool isSimplexDone(Matrix);
        unsigned getPivotRow(Matrix);
        Matrix getBasisIndices(Matrix);
    public:
        Matrix extraVariablesMatrix();
        std::vector<Matrix> initialSimplexTableau();
        LpProblem(ProblemType, Matrix, Matrix, Matrix, std::vector<restrictionType>);
        void displaySimplexTableau();
        bool isRestrictionSatisfied(Matrix, Matrix, double, restrictionType);
        bool isSolutionAdmissible(Matrix);
        void solveSimplex();
        void displayProblem();
};

#endif
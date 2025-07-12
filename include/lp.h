#ifndef LP_H
#define LP_H

#include "matrix.h"
#include <vector>
#include <string>
#include <tuple>

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
        unsigned getPivotRow(std::vector<double>, std::vector<double>, Matrix);
        Matrix getBasisIndices(Matrix);
        std::vector<std::pair<size_t, size_t>> getMaxWidth(std::vector<std::string>, Matrix, std::string mode = "matrix");
        std::vector<std::string> basisHeaders(Matrix, Matrix);
    public:
        LpProblem(void) = default;
        Matrix extraVariablesMatrix();
        std::vector<Matrix> initialSimplexTableau();
        LpProblem(ProblemType, Matrix, Matrix, Matrix, std::vector<restrictionType>);
        std::vector<std::vector<int>> getRestrictionsIndexes(Matrix);
        void displaySimplexTableau(Matrix, Matrix, Matrix, Matrix, Matrix, Matrix, Matrix);
        bool isRestrictionSatisfied(Matrix, Matrix, double, restrictionType);
        bool isSolutionAdmissible(Matrix);
        Matrix solveSimplex();
        void displayProblem();
        void addRestriction(Matrix, restrictionType, double);
        bool isProblemFeasible();
        bool isProblemBounded();
        Matrix getOptimalSolution();
        void setOptimalSolution(Matrix);
        ProblemType getType();
        Matrix getObjectiveFunction();


        // Auxiliary methods for fixing variable values
        std::tuple<LpProblem, std::vector<std::pair<unsigned, unsigned>>, double> simplifyProblem(unsigned, double);
        Matrix getSimplifiedProblemSolution(Matrix, std::vector<std::pair<unsigned, unsigned>>, unsigned, double);

        Matrix getConstraintsLHS();
        std::vector<restrictionType> getConstraintsTypes();
        Matrix getConstraintsRHS();
};

#endif
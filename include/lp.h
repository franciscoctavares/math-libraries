#ifndef LP_H
#define LP_H

#include "matrix.h"
#include "constraint.h"

#include <vector>
#include <string>
#include <tuple>

#include <limits>

#define M 1000000

enum ProblemType {
    MIN,
    MAX
};

class LpProblem {
    private:
        ProblemType type;
        Matrix objectiveFunction;
        std::vector<Constraint> constraints;

        bool isProblemCorrectlyFormulated(void);

        bool isSimplexDone(Matrix);
        unsigned getPivotRow(std::vector<double>, std::vector<double>, Matrix);
        Matrix getBasisIndices(Matrix);
        std::vector<std::pair<size_t, size_t>> getMaxWidth(std::vector<std::string>, Matrix, std::string mode = "matrix");
        std::vector<std::string> basisHeaders(Matrix, Matrix);
    public:
        LpProblem(void) = default;
        LpProblem(ProblemType, std::vector<double>, std::vector<Constraint>);

        bool isRestrictionSatisfied(std::vector<double>, int);
        bool isSolutionAdmissible(Matrix);

        Matrix extraVariablesMatrix();
        std::vector<Matrix> initialSimplexTableau();
        std::vector<std::vector<int>> getRestrictionsIndexes(Matrix);
        void displaySimplexTableau(Matrix, Matrix, Matrix, Matrix, Matrix, Matrix, Matrix);
        Matrix solveSimplex();
        void displayProblem();
        void addRestriction(Matrix, ConstraintType, double);
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
        std::vector<ConstraintType> getConstraintsTypes();
        Matrix getConstraintsRHS();
};

#endif
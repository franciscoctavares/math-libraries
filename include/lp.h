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
        Matrix optimalSolution;

        bool isProblemCorrectlyFormulated(void);

        /**
         * @brief Checks the cj - zj row see if any of the elements is positive, to see if more simplex iterations are necessary
         * 
         * @return true if all values in the cj - zj row are less than or equal to zero, false otherwise
         */
        bool isSimplexDone(Matrix);

        /**
         * @brief Given the pivot column elements(), the b column elements, and the ratios column matrix, returns the index of the pivot row
         * 
         * @return unsigned 
         */
        unsigned getPivotRow(std::vector<double>, std::vector<double>, Matrix);

        /**
         * @brief Given the extraCj row matrix(cj row minus the objective function's coefficients), return the basic variables' indexes
         * 
         * @return Matrix 
         */
        Matrix getBasisIndexes(Matrix extraCj);
        std::vector<std::pair<size_t, size_t>> getMaxWidth(std::vector<std::string>, Matrix, std::string mode = "matrix");
        std::vector<std::string> basisHeaders(Matrix, Matrix);

        Matrix getConstraintsLHS();
        std::vector<ConstraintType> getConstraintsTypes();
        Matrix getConstraintsRHS();
    public:
        /**
         * @brief Default constructor
         * 
         */
        LpProblem(void) = default;
        
        /**
         * @brief Constructs a new LP model
         * 
         * @param modelType the type of optimization problem: maximization(MAX) or minimization(MIN)
         * @param newObjectiveFunction the objective function's coefficients
         * @param newConstraints the constraints of the model
         */
        LpProblem(ProblemType modelType, std::vector<double> newObjectiveFunction, std::vector<Constraint> newConstraints);

        /**
         * @brief Checks if the provided candidate solution respects the constraint whose index is 
         *        provided by the constraintIndex argument
         * 
         * @param potentialSolution - row matrix of the candidate solution(x1, x2, x3, ... values)
         * @param constraintIndex - the index of the constraint to check
         * @return true if the constraint is satisfied and false if otherwise
         */
        bool isConstraintSatisfied(Matrix potentialSolution, int constraintIndex);

        /**
         * @brief Checks if the candidate solution satisfies all the constraints of the model
         * 
         * @param potentialSolution - the candidate solution
         * @return true if the candidate solution is feasible(satisfies all the constraints) and false otherwise
         */
        bool isSolutionAdmissible(Matrix potentialSolution);

        /**
         * @brief Builds and returns the extra variables(surplus, slack and artificial) matrix, which is used to stack horizontally to the constraints' LHS matrix
         * 
         */
        Matrix extraVariablesMatrix();

        /**
         * @brief Builds the matrices of the initial simplex tableau
         * 
         */
        std::vector<Matrix> initialSimplexTableau();

        /**
         * @brief Given the extraCj row matrix(cj row minus the objective function's coefficients), returns the extra variables' indexes in the extraCj matrix.
         *        This method is auxiliary o some other methods
         * 
         * @return std::vector<std::vector<int>> 
         */
        std::vector<std::pair<int, int>> getConstraintsIndexes(Matrix extraCj);

        /**
         * @brief Displays the current simplex tableau
         * 
         */
        void displaySimplexTableau(Matrix, Matrix, Matrix, Matrix, Matrix, Matrix, Matrix);

        /**
         * @brief Solves the LP model using the simplex method
         * 
         * @return Matrix - optimal solution
         */
        Matrix solveSimplex();

        /**
         * @brief Displays the optimal solution to the model
         * 
         */
        void displayProblem(Matrix optimalSolution);
        
        /**
         * @brief Adds a new constraint to the LP model
         * 
         * @param newConstraint - the new constraint to be added
         */
        void addConstraint(Constraint newConstraint);
        
        /**
         * @brief Checks if the output of the simplex solver is a feasible solution, that is, if all artificial variables(if present)
         *        were eliminated
         * 
         * @param candidateSolution - the potential feasible solution 
         * @return true if the solution is feasible, false otherwise
         */
        bool isProblemFeasible(Matrix candidateSolution);

        /**
         * @brief Checks if the problem is unbounded, using the output of the simplex solver.
         * 
         * @param candidateSolution - the output of the simplex solver
         * @return true if the problem is unbounded, false otherwise
         */
        bool isProblemBounded(Matrix candidateSolution);


        Matrix getOptimalSolution();
        void setOptimalSolution(Matrix);
        ProblemType getType();
        Matrix getObjectiveFunction();


        // Auxiliary methods for fixing variable values
        std::tuple<LpProblem, std::vector<std::pair<unsigned, unsigned>>, double> simplifyProblem(unsigned, double);
        Matrix getSimplifiedProblemSolution(Matrix, std::vector<std::pair<unsigned, unsigned>>, unsigned, double);
};

#endif
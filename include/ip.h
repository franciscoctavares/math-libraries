#ifndef IP_H
#define IP_H

#include "lp.h"
#include <vector>

enum nodeStatus {
    EVALUATING,
    FATHOMED,
};

struct BBNode {
    LpProblem problem;
    BBNode* child1;
    BBNode* child2;
    nodeStatus status;
};

class BranchAndBoundNode {
    private:
        LpProblem problem;
        BranchAndBoundNode* child1;
        BranchAndBoundNode* child2;
        nodeStatus status;
    public:
        BranchAndBoundNode(void) = default;
        BranchAndBoundNode(LpProblem);
        Matrix solveCurrentNode();
        /**
         * @brief Determines if the presented solution to an Integer Programming Problem is composed entirely of integer variables
         * 
         * @return true, if the solution is composed only of integers, false otherwise.
         */
        bool isSolutionWhole(Matrix);
        unsigned getBranchVariableIndex(Matrix);
        Matrix solveChildrenNodes(unsigned);
        Matrix buildNewRestriction(unsigned, restrictionType, double);
        Matrix compareChildrenSolutions(Matrix, Matrix, Matrix);
        std::pair<int, double> canProblemBeSimplified(Matrix, restrictionType, double);
};

class BranchAndBoundTree {
    private:
        BranchAndBoundNode headNode;
    public:
        BranchAndBoundTree(LpProblem);
        BranchAndBoundNode getHeadNode();

};

class IpProblem {
    private:
        LpProblem startingProblem;
        Matrix optimalSolution;
    public:
        IpProblem(LpProblem);
        void solveBranchAndBound();
};

#endif
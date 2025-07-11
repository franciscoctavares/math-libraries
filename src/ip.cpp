#include "../include/ip.h"

#include <cmath>
#include <iostream>

BranchAndBoundNode::BranchAndBoundNode(LpProblem nodeProblem) {
    problem = nodeProblem;
    child1 = nullptr;
    child2 = nullptr;
    status = EVALUATING;
}

Matrix BranchAndBoundNode::solveCurrentNode() {
    Matrix solution = problem.solveSimplex();
    problem.displayProblem();
    return solution;
}

bool BranchAndBoundNode::isSolutionWhole(Matrix solution) {
    Matrix currentSolution = problem.getOptimalSolution();
    for(int i = 0; i < currentSolution.columns(); i++) {
        if(currentSolution.getElement(0, i) != floor(currentSolution.getElement(0, i))) return false;
    }
    return true;
}

unsigned BranchAndBoundNode::getBranchVariableIndex(Matrix solution) {
    for(unsigned i = 0; i < solution.columns(); i++) {
        if(solution.getElement(0, i) != floor(solution.getElement(0, i))) return i;
    }
    return 0;
}

Matrix BranchAndBoundNode::solveChildrenNodes(unsigned n) {
    //std::cout << "This is problem " << n << std::endl;
    //problem.displaySimplexTableau();
    Matrix currentSol = problem.solveSimplex();
    problem.displayProblem();

    if(isSolutionWhole(currentSol)) {
        std::cout << "The solution to the problem " << n << " is whole." << std::endl;
        return currentSol;
    }
    else if(!problem.isProblemBounded() || !problem.isProblemFeasible()) return currentSol;

    // calculate branching variable and its bounds for left and right child
    unsigned branchIndex = getBranchVariableIndex(currentSol);
    double value = currentSol.getElement(0, branchIndex);
    double valueChild1 = floor(value);
    double valueChild2 = ceil(value);
    //std::cout << "The branching variable is x" << branchIndex + 1 << ". ";

    // Left child
    LpProblem child1Aux = problem;
    Matrix newLhs = basisVector(currentSol.columns(), branchIndex).transpose();
    child1Aux.addRestriction(newLhs, LESS_THAN_OR_EQUAL, valueChild1);
    child1 = new BranchAndBoundNode(child1Aux);
    Matrix child1Sol = child1->solveChildrenNodes(n + 1);

    // Right child
    LpProblem child2Aux = problem;
    newLhs = basisVector(currentSol.columns(), branchIndex).transpose();
    child2Aux.addRestriction(newLhs, GREATER_THAN_OR_EQUAL, valueChild2);
    //std::cout << "This is child2" << std::endl;
    child2 = new BranchAndBoundNode(child2Aux);
    Matrix child2Sol = child2->solveChildrenNodes(n + 2);
    //delete child2;

    // Calculate the better solution
    Matrix betterSol = compareChildrenSolutions(child1Sol, child2Sol, problem.getObjectiveFunction());
    //betterSol.displayMatrix();

    delete child1;
    delete child2;
    return betterSol;
}

Matrix BranchAndBoundNode::buildNewRestriction(unsigned varIndex, restrictionType restType, double rhs) {
    Matrix lhs = basisVector(problem.getOptimalSolution().columns(), varIndex).transpose();

}

Matrix BranchAndBoundNode::compareChildrenSolutions(Matrix child1Sol, Matrix child2Sol, Matrix objectiveFunction) {
    bool child1Valid, child2Valid;
    if(child1Sol.rows() == 1 && child1Sol.columns() == 1 && child1Sol.getElement(0, 0) == 0) child1Valid = false;
    else if(child1Sol.rows() == 1 && child1Sol.columns() == 1 && child1Sol.getElement(0, 0) == INFINITY) child1Valid = false;
    else child1Valid = true;

    if(child2Sol.rows() == 1 && child2Sol.columns() == 1 && child2Sol.getElement(0, 0) == 0) child2Valid = false;
    else if(child2Sol.rows() == 1 && child2Sol.columns() == 1 && child2Sol.getElement(0, 0) == INFINITY) child2Valid = false;
    else child2Valid = true;

    if(child1Valid && !child2Valid) return child1Sol;
    else if(!child1Valid && child2Valid) return child2Sol;
    else if(!child1Valid && !child2Valid) return child1Sol;
    else {
        std::cout << "Both solutions are valid" << std::endl;
        double child1Result = child1Sol.dotProduct(objectiveFunction);
        double child2Result = child2Sol.dotProduct(objectiveFunction);
        if(problem.getType() == MAX) return (child1Result > child2Result) ? child1Sol : child2Sol;
        else if(problem.getType() == MIN) return (child1Result < child2Result) ? child1Sol : child2Sol;
    }
}

BranchAndBoundTree::BranchAndBoundTree(LpProblem startingProblem) {
    headNode = BranchAndBoundNode(startingProblem);
}

BranchAndBoundNode BranchAndBoundTree::getHeadNode() {
    return headNode;
}

IpProblem::IpProblem(LpProblem problem) {
    startingProblem = problem;
}

void IpProblem::solveBranchAndBound() {
    BranchAndBoundTree tree(startingProblem);
    //Matrix solution = tree.getHeadNode().solveCurrentNode();
    //unsigned lol = tree.getHeadNode().getBranchVariableIndex(solution);
    Matrix bestSolution = tree.getHeadNode().solveChildrenNodes(0);
    bestSolution.displayMatrix();
}
#include "../include/lp.h"

LpProblem::LpProblem(Matrix objCoeffs, Matrix restLHS, Matrix restRHS, std::vector<restrictionType> restType) {
    objectiveFunction = objCoeffs; // assumes objectiveFunction is a row matrix
    restrictionsLHS = restLHS;     // assumes restrictionsLHS is a matrix
    restrictionsRHS = restRHS;     // assumes restrictionsRHS is a column matrix
    restrictionsTypes = restType;
}

bool LpProblem::isRestrictionSatisfied(Matrix potentialSolution, Matrix restLHS, double restRHS, restrictionType restType) {
    double value = potentialSolution.dotProduct(restLHS);
    if(restType == LESS_THAN_OR_EQUAL) return (value <= restRHS) ? true : false;
    else if(restType == EQUAL) return (value == restRHS) ? true : false;
    else if(restType == GREATER_THAN_OR_EQUAL) return (value >= restRHS) ? true : false;
}

bool LpProblem::isSolutionAdmissible(Matrix potentialSolution) {
    for(int i = 0; i < restrictionsLHS.rows(); i++) {
        if(!isRestrictionSatisfied(potentialSolution, restrictionsLHS.getRow(i), restrictionsRHS.getElement(i, 0), restrictionsTypes[i])) return false;
    }
    return true;
}
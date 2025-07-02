#include "../include/lp.h"
#include <cmath>
#include <iostream>
#include <iomanip>

LpProblem::LpProblem(ProblemType probType, Matrix objCoeffs, Matrix restLHS, Matrix restRHS, std::vector<restrictionType> restType) {
    type = probType,
    objectiveFunction = objCoeffs; // assumes objectiveFunction is a row matrix
    restrictionsLHS = restLHS;
    restrictionsRHS = restRHS;     // assumes restrictionsRHS is a column matrix
    restrictionsTypes = restType;
}

void LpProblem::displaySimplexTableau() {
    // objective function coefficients
    std::cout << "|";
    for(int i = 0; i < objectiveFunction.columns(); i++) {
        if(objectiveFunction.getElement(0, i) < 0) std::cout << "-";
        else std::cout << "+";
        std::cout << std::setprecision(3) << std::fixed << fabs(objectiveFunction.getElement(0, i));
        if(i < objectiveFunction.columns() - 1) std::cout << " ";
    }
    std::cout << "|" << std::endl << std::endl;

    // restrictions
    for(int i = 0; i < restrictionsLHS.rows(); i++) {
        std::cout << "|";
        for(int j = 0; j < restrictionsLHS.columns(); j++) {
            if(restrictionsLHS.getElement(i, j) < 0) std::cout << "-";
            else std::cout << "+";
            std::cout << std::setprecision(3) << std::fixed << fabs(restrictionsLHS.getElement(i, j));
            if(j < restrictionsLHS.columns() - 1) std::cout << " ";
        }
        std::cout << "|";

        if(restrictionsTypes[i] == LESS_THAN_OR_EQUAL) std::cout << " <= ";
        else if(restrictionsTypes[i] == EQUAL) std::cout << " = ";
        else if(restrictionsTypes[i] == GREATER_THAN_OR_EQUAL) std::cout << " >= ";

        std::cout << "|";

        if(restrictionsRHS.getElement(i, 0) < 0) std::cout << "-";
        else std::cout << "+";

        std::cout << std::setprecision(3) << std::fixed << fabs(restrictionsRHS.getElement(i, 0));
        std::cout << "|";

        std::cout << std::endl;
    }
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

    // non negativity
    for(int i = 0; i < potentialSolution.columns(); i++) {
        if(potentialSolution.getElement(0, i) < 0) return false;
    }

    return true;
}

std::vector<Matrix> LpProblem::initialSimplexTableau() {
    // to determine how many extra columns the simplex tableau will need
    unsigned slack_surplus_variables = 0;
    for(int i = 0; i < restrictionsLHS.rows(); i++) {
        if(restrictionsTypes[i] == LESS_THAN_OR_EQUAL) slack_surplus_variables += 1;
        else if(restrictionsTypes[i] == EQUAL) slack_surplus_variables += 0;
        else if(restrictionsTypes[i] == GREATER_THAN_OR_EQUAL) slack_surplus_variables += 1;
    }

    unsigned artificial_variables = 0;
    for(int i = 0; i < restrictionsLHS.rows(); i++) {
        if(restrictionsTypes[i] == GREATER_THAN_OR_EQUAL) artificial_variables += 1;
    }

    // extra variables matrix
    std::vector<Matrix> extraStuff;
    for(int i = 0; i < restrictionsLHS.rows(); i++) {
        if(restrictionsTypes[i] == LESS_THAN_OR_EQUAL) extraStuff.push_back(basisVector(restrictionsLHS.rows(), i));
        else if(restrictionsTypes[i] == EQUAL) continue;
        else if(restrictionsTypes[i] == GREATER_THAN_OR_EQUAL) extraStuff.push_back(basisVector(restrictionsLHS.rows(), i) * -1);
    }
    // artificial variables
    for(int i = 0; i < restrictionsLHS.rows(); i++) {
        if(restrictionsTypes[i] == GREATER_THAN_OR_EQUAL) extraStuff.push_back(basisVector(restrictionsLHS.rows(), i));
    }

    Matrix extraMatrix = extraStuff[0];
    for(int i = 1; i < extraStuff.size(); i++) {
        extraMatrix = extraMatrix.stackHorizontal(extraStuff[i]);
    }

    Matrix simplexTableau = restrictionsLHS.stackHorizontal(extraMatrix);

    //simplexTableau.displayMatrix();
    //std::cout << std::endl;

    Matrix b = restrictionsRHS;
    //b.displayMatrix();
    //std::cout << std::endl;

    Matrix cj = objectiveFunction;
    if(type == MIN) cj = cj * -1;
    std::vector<Matrix> extraAux;
    std::vector<double> zero = {0};
    Matrix mZero(zero, 1, 1);
    std::vector<double> bigM = {M};
    Matrix mBigM(bigM, 1, 1);
    for(int i = 0; i < restrictionsLHS.rows(); i++) {
        if(restrictionsTypes[i] == LESS_THAN_OR_EQUAL) extraAux.push_back(mZero);
        else if(restrictionsTypes[i] == EQUAL) extraAux.push_back(mBigM);
        else if(restrictionsTypes[i] == GREATER_THAN_OR_EQUAL) extraAux.push_back(mBigM);
        //cj.displayMatrix();
        //std::cout << cj.columns() << std::endl;
    }
    for(int i = 0; i < extraAux.size(); i++) {
        cj = cj.stackHorizontal(extraAux[i]);
    }
    //cj.displayMatrix();
    //std::cout << std::endl;

    Matrix zj = zeros(1, simplexTableau.columns());
    Matrix cj_minus_zj = zeros(1, simplexTableau.columns());

    unsigned n_variables = objectiveFunction.columns();
    unsigned extra_variables = cj.columns() - n_variables;

    //Matrix basisVariablesIndices = zeros(cb.rows(), 1);

    std::vector<double> basisThing;
    //Matrix basisVariablesIndices = zeros(cb.rows(), 1);
    std::vector<double> basisIndices;
    for(int i = 0; i < restrictionsLHS.rows(); i++) {
        if(restrictionsTypes[i] == LESS_THAN_OR_EQUAL) {
            basisThing.push_back(0.0);
        }
        else if(restrictionsTypes[i] == EQUAL) {
            basisThing.push_back(-1 * M);
        }
        else if(restrictionsTypes[i] == GREATER_THAN_OR_EQUAL) {
            basisThing.push_back(-1 * M);
        }
        basisIndices.push_back(n_variables + i);
    }
    Matrix cb(basisThing, basisThing.size(), 1);
    //cb.displayMatrix();
    //std::cout << std::endl;
    Matrix basisVariablesIndices(basisIndices, cb.rows(), 1);

    std::vector<Matrix> outputMatrices = {simplexTableau, b, cj, basisVariablesIndices, cb};
    return outputMatrices;
}

bool LpProblem::isSimplexDone(Matrix cj_minus_zj) {
    for(int i = 0; i < cj_minus_zj.columns(); i++) {
        if(cj_minus_zj.getElement(0, i) > 0) return false;
    }
    return true;
}

unsigned LpProblem::getPivotRow(Matrix ratios) {
    double minValue = M;
    unsigned minIndex = 0;
    for(int i = 0; i < ratios.rows(); i++) {
        if(ratios.getElement(i, 0) < minValue && ratios.getElement(i, 0) >= 0) {
            minValue = ratios.getElement(i, 0);
            minIndex = i;
        }
    }
    return minIndex;
}

Matrix LpProblem::solveSimplex() {
    Matrix pivots = zeros(1, 2);
    std::vector<Matrix> things = initialSimplexTableau();
    Matrix simplexTableau = things[0];
    Matrix b = things[1];
    Matrix cj = things[2];
    Matrix basisIndices = things[3];
    Matrix cb = things[4];

    Matrix zj = zeros(1, simplexTableau.columns());
    Matrix cj_minus_zj = zeros(1, simplexTableau.columns());

    for(int i = 0; i < simplexTableau.columns(); i++) {
        zj.setElement(0, i, cb.dotProduct(simplexTableau.getColumn(i)));
    }

    //zj.displayMatrix();
    //std::cout << std::endl;

    cj_minus_zj = cj - zj;
    //cj_minus_zj.displayMatrix();
    //std::cout << std::endl;

    //cj.displayMatrix();
    //std::cout << std::endl;

    //b.displayMatrix();
    //std::cout << std::endl;


    while(!isSimplexDone(cj_minus_zj)) {
        pivots.setElement(0, 1, cj_minus_zj.maxValueIndex());
        Matrix ratios = zeros(restrictionsLHS.rows(), 1);
        ratios = b.pointDivision(simplexTableau.getColumn(pivots.getElement(0, 1)));
        //ratios.displayMatrix();
        //std::cout << std::endl;

        pivots.setElement(0, 0, getPivotRow(ratios));
        //pivots.displayMatrix();
        //std::cout << std::endl;

        unsigned oldBasis = pivots.getElement(0, 0);
        unsigned newBasis = pivots.getElement(0, 1);
        //pivots.displayMatrix();
        //std::cout << std::endl;
        basisIndices.setElement(oldBasis, 0, newBasis);
        cb.setElement(oldBasis, 0, cj.getElement(0, newBasis));
        //cb.displayMatrix();
        //std::cout << std::endl;

        Matrix newRow = simplexTableau.getRow(oldBasis) * (1 / simplexTableau.getElement(oldBasis, newBasis));
        b.setElement(oldBasis, 0, b.getElement(oldBasis, 0) / simplexTableau.getElement(oldBasis, newBasis));
        
        //b.displayMatrix();
        //std::cout << std::endl;

        //std::cout << b.getElement(oldBasis, 0) << " " << simplexTableau.getElement(oldBasis, newBasis) << std::endl;
        //newRow.displayMatrix();
        //std::cout << std::endl;

        simplexTableau = simplexTableau.setRow(oldBasis, newRow);
        //simplexTableau.displayMatrix();
        //std::cout << std::endl;

        for(int i = 0; i < simplexTableau.rows(); i++) {
            if(i == oldBasis) continue;
            else {
                double factor = simplexTableau.getElement(i, newBasis);
                simplexTableau = simplexTableau.rowOperation(oldBasis, i, -1 * factor);
                b = b.rowOperation(oldBasis, i, -1 * factor);
                //b.setElement(i, 0, factor * b.getElement(oldBasis, 0));
            }
        }

        //simplexTableau.displayMatrix();
        //std::cout << std::endl;

        //b.displayMatrix();
        //std::cout << std::endl;

        for(int i = 0; i < simplexTableau.columns(); i++) {
            zj.setElement(0, i, cb.dotProduct(simplexTableau.getColumn(i)));
        }
        cj_minus_zj = cj - zj;
        //cj_minus_zj.displayMatrix();
        //std::cout << std::endl;

    }
    //std::cout << "Pivot column is " << unsigned(pivots.getElement(0, 1)) << std::endl;

    //Matrix ratios = zeros(restrictionsLHS.rows(), 1);
    //ratios = b.pointDivision(simplexTableau.getColumn(pivots.getElement(0, 1)));
    //ratios.displayMatrix();
    //std::cout << std::endl;

    //pivots.setElement(0, 0, ratios.minValueIndex());
    //std::cout << "Pivot row is " << unsigned(pivots.getElement(0, 0)) << std::endl;

    //auxiliar

    Matrix solution = zeros(1, objectiveFunction.columns());

    //simplexTableau.displayMatrix();
    //std::cout << std::endl;

    //b.displayMatrix();
    //std::cout << std::endl;

    for(int k = 0; k < basisIndices.rows(); k++) {
        solution.setElement(0, basisIndices.getElement(k, 0), b.getElement(k, 0));
    }

    optimalSolution = solution;

    return solution;
}

void LpProblem::displayProblem() {
    std::cout << "The optimal solution is (";
    for(int i = 0; i < optimalSolution.columns(); i++) {
        std::cout << "x" << i + 1;
        if(i < optimalSolution.columns() - 1) std::cout << ", ";
    }
    std::cout << ") = (";
    for(int i = 0; i < optimalSolution.columns(); i++) {
        std::cout << std::setprecision(3) << std::fixed << optimalSolution.getElement(0, i);
        if(i < optimalSolution.columns() - 1) std::cout << ", ";
    }
    std::cout << ")" << std::endl;
}
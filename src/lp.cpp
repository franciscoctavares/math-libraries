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
    if(restType == LESS_THAN_OR_EQUAL) return (value <= restRHS) ? true : false;    // <=
    else if(restType == EQUAL) return (value == restRHS) ? true : false;            // =
    else return (value >= restRHS) ? true : false;                                  // >=
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

Matrix LpProblem::getBasisIndices(Matrix extraCj) {
    std::vector<std::vector<int>> restrictionsIndices;
    for(int i = 0; i < restrictionsTypes.size(); i++) restrictionsIndices.push_back({-1, -1});


    unsigned n_slack_surplus_variables = 0;
    unsigned totalExtraVariables = extraCj.columns();
    unsigned currentCoefficient = 0;
    for(int i = 0; i < restrictionsTypes.size(); i++) {
        if(restrictionsTypes[i] == LESS_THAN_OR_EQUAL) {
            restrictionsIndices[i][0] = currentCoefficient;
            currentCoefficient++;
            restrictionsIndices[i][1] = -2; // nas restrições <= não há variáveis artificiais
            n_slack_surplus_variables++;
        }
        else if(restrictionsTypes[i] == GREATER_THAN_OR_EQUAL) {
            restrictionsIndices[i][0] = currentCoefficient;
            currentCoefficient++;
            n_slack_surplus_variables++;
        }
        else if(restrictionsTypes[i] == EQUAL) {
            restrictionsIndices[i][0] = -2; // nas restrições = não há coeficiente do zero
        }
    }

    std::vector<unsigned> artificialRestrictions;
    for(int i = 0; i < restrictionsTypes.size(); i++) {
        if(restrictionsIndices[i][1] == -2) continue;
        else if(restrictionsIndices[i][1] == -1) artificialRestrictions.push_back(i);
    }

    for(int i = n_slack_surplus_variables; i < totalExtraVariables; i++) {
        if(artificialRestrictions.size() > 0) {
            restrictionsIndices[artificialRestrictions[0]][1] = i;
            artificialRestrictions.erase(artificialRestrictions.begin());
        }
    }

    for(int i = 0; i < restrictionsIndices.size(); i++) {
        if(restrictionsIndices[i][0] != -2) restrictionsIndices[i][0] += objectiveFunction.columns();
        if(restrictionsIndices[i][1] != -2) restrictionsIndices[i][1] += objectiveFunction.columns();
        //std::cout << "[" << restrictionsIndices[i][0] << ", " << restrictionsIndices[i][1] << "]" << std::endl;
    }

    std::vector<double> basisIndices;
    for(int i = 0; i < restrictionsTypes.size(); i++) {
        if(restrictionsTypes[i] == LESS_THAN_OR_EQUAL) basisIndices.push_back(restrictionsIndices[i][0]);
        else if(restrictionsTypes[i] == GREATER_THAN_OR_EQUAL) basisIndices.push_back(restrictionsIndices[i][1]);
        else if(restrictionsTypes[i] == EQUAL) basisIndices.push_back(restrictionsIndices[i][1]);
    }

    return Matrix(basisIndices, restrictionsTypes.size(), 1);
}

std::vector<Matrix> LpProblem::initialSimplexTableau() {
    Matrix simplexTableau = restrictionsLHS;
    simplexTableau.stackHorizontal(extraVariablesMatrix());
    Matrix b = restrictionsRHS;
    Matrix cj = objectiveFunction;
    
    if(type == MIN) cj = cj * -1;

    std::vector<double> aux;
    for(int i = 0; i < restrictionsLHS.rows(); i++) {
        if(restrictionsTypes[i] == LESS_THAN_OR_EQUAL || restrictionsTypes[i] == GREATER_THAN_OR_EQUAL) aux.push_back(0.0);
    }
    for(int i = 0; i < restrictionsLHS.rows(); i++) {
        if(restrictionsTypes[i] == EQUAL || restrictionsTypes[i] == GREATER_THAN_OR_EQUAL) {
            if(type == MAX) aux.push_back(-1 * M);
            else if(type == MIN) aux.push_back(M);
        }
    }

    Matrix extraCj(aux, 1, aux.size());

    Matrix basisIndicesAux = getBasisIndices(extraCj);
    cj.stackHorizontal(extraCj);

    std::vector<double> basisThing;
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
    }
    Matrix cb(basisThing, basisThing.size(), 1);

    return {simplexTableau, b, cj, basisIndicesAux, cb};
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
        if(ratios.getElement(i, 0) < minValue && ratios.getElement(i, 0) >= 0 && ratios.getElement(i, 0) != M) {
            minValue = ratios.getElement(i, 0);
            minIndex = i;
        }
    }
    return minIndex;
}

Matrix LpProblem::extraVariablesMatrix() {
    unsigned nVariables = 0;
    unsigned slack_surplus_variables = 0;
    unsigned artificial_variables = 0;

    std::vector<std::vector<double>> pairs;

    for(int i = 0; i < restrictionsTypes.size(); i++) {
        if(restrictionsTypes[i] == LESS_THAN_OR_EQUAL) {
            nVariables += 1;
            slack_surplus_variables += 1;
        }
        else if(restrictionsTypes[i] == EQUAL) {
            nVariables += 1;
            artificial_variables += 1;
        }
        else if(restrictionsTypes[i] == GREATER_THAN_OR_EQUAL) {
            nVariables += 2;
            slack_surplus_variables += 1;
            artificial_variables += 1;
        }
    }

    for(int i = 0; i < restrictionsTypes.size(); i++) {
        if(restrictionsTypes[i] == LESS_THAN_OR_EQUAL) {
            pairs.push_back({double(i), 1.0});
        }
        else if(restrictionsTypes[i] == EQUAL) {
            continue;
        }
        else if(restrictionsTypes[i] == GREATER_THAN_OR_EQUAL) {
            pairs.push_back({double(i), -1.0});
        }
    }

    for(int i = 0; i < restrictionsTypes.size(); i++) {
        if(restrictionsTypes[i] == EQUAL || restrictionsTypes[i] == GREATER_THAN_OR_EQUAL) {
            //pairs.push_back({double(slack_surplus_variables - 1 + i), 1.0});
            pairs.push_back({double(i), 1.0});
        }
        else continue;
    }

    Matrix aux = zeros(restrictionsLHS.rows(), nVariables);

    for(int i = 0; i < pairs.size(); i++) {
        aux = aux.setColumn(i, basisVector(restrictionsLHS.rows(), pairs[i][0]) * pairs[i][1]);
    }

    return aux;
}

void LpProblem::solveSimplex() {
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
    cj_minus_zj = cj - zj;

    while(!isSimplexDone(cj_minus_zj)) {
        pivots.setElement(0, 1, cj_minus_zj.maxValueIndex());
        Matrix ratios = zeros(restrictionsLHS.rows(), 1);
        ratios = b.pointDivision(simplexTableau.getColumn(pivots.getElement(0, 1)));

        pivots.setElement(0, 0, getPivotRow(ratios));

        unsigned oldBasis = pivots.getElement(0, 0);
        unsigned newBasis = pivots.getElement(0, 1);

        if(cb.getElement(oldBasis, 0) == M || cb.getElement(oldBasis, 0) == -1 * M) {
            unsigned artificial_index = basisIndices.getElement(oldBasis, 0);

            simplexTableau = simplexTableau.removeColumn(artificial_index);
            cj = cj.removeColumn(artificial_index);
            zj = zj.removeColumn(artificial_index);
            cj_minus_zj = cj_minus_zj.removeColumn(artificial_index);

            for(int i = 0; i < basisIndices.rows(); i++) {
                if(basisIndices.getElement(i, 0) > artificial_index) basisIndices.setElement(i, 0, basisIndices.getElement(i, 0) - 1.0);
            }
        }

        basisIndices.setElement(oldBasis, 0, newBasis);
        cb.setElement(oldBasis, 0, cj.getElement(0, newBasis));

        Matrix newRow = simplexTableau.getRow(oldBasis) * (1 / simplexTableau.getElement(oldBasis, newBasis));
        b.setElement(oldBasis, 0, b.getElement(oldBasis, 0) / simplexTableau.getElement(oldBasis, newBasis));

        simplexTableau = simplexTableau.setRow(oldBasis, newRow);

        for(int i = 0; i < simplexTableau.rows(); i++) {
            if(i == oldBasis) continue;
            else {
                double factor = simplexTableau.getElement(i, newBasis);
                simplexTableau.rowOperation(oldBasis, i, -1 * factor);
                b.rowOperation(oldBasis, i, -1 * factor);
            }
        }

        for(int i = 0; i < simplexTableau.columns(); i++) {
            zj.setElement(0, i, cb.dotProduct(simplexTableau.getColumn(i)));
        }
        cj_minus_zj = cj - zj;

    }

    Matrix solution = zeros(1, objectiveFunction.columns());

    for(int k = 0; k < basisIndices.rows(); k++) {
        solution.setElement(0, basisIndices.getElement(k, 0), b.getElement(k, 0));
    }

    optimalSolution = solution;
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
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

std::vector<std::vector<int>> LpProblem::getRestrictionsIndexes(Matrix extraCj) {
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
    return restrictionsIndices;
}

size_t LpProblem::getMaxWidth(std::vector<std::string> headers, Matrix matrix) {
    size_t maxWidth = 0;
    for(double val: matrix.getElements()) {
        std::string str = std::to_string(val);
        // Trim trailing zeroes for nicer formatting (optional)
        str.erase(str.find_last_not_of('0') + 1, std::string::npos);
        if (str.back() == '.') str.pop_back(); // remove trailing dot if needed
        maxWidth = std::max(maxWidth, str.length());
    }
    return maxWidth;
}

void LpProblem::displaySimplexTableau(Matrix tableau, Matrix cb, Matrix basisIndexes, Matrix cj, Matrix b, Matrix zj, Matrix cj_minus_zj) {
    Matrix extraCj = cj.subMatrix(0, 0, objectiveFunction.columns(), cj.columns() - 1);
    std::vector<std::vector<int>> restrictionsIndexes = getRestrictionsIndexes(extraCj);

    // compute the maximum width of all numbers to be printed
    size_t maxWidth = 0;


    for(double val: tableau.getElements()) {
        std::string str = std::to_string(val);
        // Trim trailing zeroes for nicer formatting (optional)
        str.erase(str.find_last_not_of('0') + 1, std::string::npos);
        if (str.back() == '.') str.pop_back(); // remove trailing dot if needed
        maxWidth = std::max(maxWidth, str.length());
    }
    for(double val: cj.getElements()) {
        std::string str = std::to_string(val);
        // Trim trailing zeroes for nicer formatting (optional)
        str.erase(str.find_last_not_of('0') + 1, std::string::npos);
        if (str.back() == '.') str.pop_back(); // remove trailing dot if needed
        maxWidth = std::max(maxWidth, str.length());
    }
    for(double val: zj.getElements()) {
        std::string str = std::to_string(val);
        // Trim trailing zeroes for nicer formatting (optional)
        str.erase(str.find_last_not_of('0') + 1, std::string::npos);
        if (str.back() == '.') str.pop_back(); // remove trailing dot if needed
        maxWidth = std::max(maxWidth, str.length());
    }
    for(double val: cj_minus_zj.getElements()) {
        std::string str = std::to_string(val);
        // Trim trailing zeroes for nicer formatting (optional)
        str.erase(str.find_last_not_of('0') + 1, std::string::npos);
        if (str.back() == '.') str.pop_back(); // remove trailing dot if needed
        maxWidth = std::max(maxWidth, str.length());
    }

    std::vector<std::string> variables;
    for(int i = 0; i < objectiveFunction.columns(); i++) {
        variables.push_back("x" + std::to_string(i + 1));
    }
    for(int k = 0; k < restrictionsIndexes.size(); k++) {
        if(restrictionsIndexes[k][0] > 0) variables.push_back("s" + std::to_string(k + 1));
    }
    for(int k = 0; k < restrictionsIndexes.size(); k++) {
        if(restrictionsIndexes[k][1] > 0) variables.push_back("a" + std::to_string(k + 1));
    }

    std::ostringstream oss;

    // linha do cj
    std::cout << "|";
    oss << "|";
    
    
    for(int i = 0; i < objectiveFunction.columns(); i++) {
        std::cout << std::setw(maxWidth) << "x" << i + 1;
        oss << std::setw(maxWidth) << "x" << i + 1;
    }
    for(int j = 0; j < restrictionsIndexes.size(); j++) {
        if(restrictionsIndexes[j][0] > 0) {
            std::cout << std::setw(maxWidth) << "s" << j + 1;
            oss << std::setw(maxWidth) << "s" << j + 1;
        }
    }
    for(int k = 0; k < restrictionsIndexes.size(); k++) {
        if(restrictionsIndexes[k][1] > 0) {
            std::cout << std::setw(maxWidth) << "a" << k + 1;
            oss << std::setw(maxWidth) << "a" << k + 1;
        }
    }
    

    std::cout << "|" << std::endl;
    std::cout << "|";
    for(int i = 0; i < cj.columns(); i++) {
        std::cout << std::setw(maxWidth + 1) << cj.getElement(0, i);
    }
    std::cout << "|";
    oss << "|";

    std::string line = oss.str();
    size_t line_length = line.length();
    std::cout << std::endl;
    for(int i = 0; i < line_length; i++) std::cout << "-";

    // basis variables and its coefficients

    std::vector<std::string> basisVariables;
    for(int i = 0; i < basisIndexes.getElements().size(); i++) {

    }

    std::cout << std::endl;
}

bool LpProblem::isRestrictionSatisfied(Matrix potentialSolution, Matrix restLHS, double restRHS, restrictionType restType) {
    double value = potentialSolution.dotProduct(restLHS);
    if(restType == LESS_THAN_OR_EQUAL) return (value <= restRHS) ? true : false;    // <=
    else if(restType == EQUAL) return (value == restRHS) ? true : false;            // =
    else return (value >= restRHS) ? true : false;                                  // >=
}

bool LpProblem::isSolutionAdmissible(Matrix potentialSolution) {
    for(int i = 0; i < restrictionsLHS.rows(); i++) {
        if(!isRestrictionSatisfied(potentialSolution, Matrix(restrictionsLHS.getRow(i), 1, restrictionsLHS.columns()), restrictionsRHS.getElement(i, 0), restrictionsTypes[i])) return false;
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
            if(type == MIN) aux.push_back(-1 * M);
            else if(type == MAX) aux.push_back(M);
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

unsigned LpProblem::getPivotRow(std::vector<double> simplexAux, std::vector<double> bAux, Matrix ratios) {
    /*
    //b.displayMatrix();
    //std::cout << std::endl;
    
    //simplex.displayMatrix();
    //std::cout << std::endl;

    ratios.displayMatrix();
    std::cout << std::endl;


    double minValue = M;
    unsigned minIndex = 0;
    
    std::vector<unsigned> validIndices;
    for(int i = 0; i < ratios.rows(); i++) validIndices.push_back(i);
    for(int i = 0; i < ratios.rows(); i++) {
        if(std::isinf(ratios.getElement(i, 0))) validIndices[i] = -1;
        else if(ratios.getElement(i, 0) < 0) validIndices[i] = -1;
    }

    std::vector<double> validValues;
    for(int i = 0; i < validIndices.size(); i++) {
        if(validIndices[i] != -1) validValues.push_back(ratios.getElement(i, 0));
        else validValues.push_back(-1.0);
    }

    for(int i = 0; i < validIndices.size(); i++) std::cout << validValues[i] << " ";
    std::cout << std::endl;

    for(int i = 0; i < validValues.size(); i++) {
        if(validValues[i] == -1.0) continue;
        if(validValues[minIndex] == -1.0) minIndex = i;
        else if(validValues[i] < validValues[minIndex]) minIndex;
    }

    std::cout << "minIndex = " << minIndex << std::endl;
    return minIndex;

    for(int i = 0; i < ratios.rows(); i++) {
        //if(simplex.getElement(i, 0) <= 0 || ratios.getElement(i, 0) == 0) continue;
        if(ratios.getElement(i, 0) < 0) continue;
        else if(std::isnan(ratios.getElement(i, 0))) continue;

        if(ratios.getElement(i, 0) < minValue && ratios.getElement(i, 0) >= 0 && ratios.getElement(i, 0) != M) {
            minValue = ratios.getElement(i, 0);
            minIndex = i;
        }
    }



    return minIndex;
    */
    //ratios.displayMatrix();
    //std::cout << std::endl;
    //for(int i = 0; i < simplexAux.size(); i++) std::cout << simplexAux[i] << " " << bAux[i] << std::endl;

    ratios.displayMatrix();
    std::cout << std::endl;

    double minValue = M;
    unsigned minIndex = 0;
    unsigned invalid = 0;
    for(int i = 0; i < ratios.rows(); i++) {
        if(simplexAux[i] <= 0) invalid++;
    }
    if(invalid == ratios.rows()) return -1;
    for(int i = 0; i < ratios.rows(); i++) {
        if(ratios.getElement(i, 0) < minValue && ratios.getElement(i, 0) > 0 && ratios.getElement(i, 0) != M) {
            minValue = ratios.getElement(i, 0);
            minIndex = i;
        }
    }
    std::cout << "minValue = " << ratios.getElement(minIndex, 0) << " and minIndex = " << minIndex << std::endl << std::endl;
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

    
    displaySimplexTableau(simplexTableau, cb, basisIndices, cj, b, zj, cj_minus_zj);
    return zeros(1, objectiveFunction.columns());

    unsigned n_surplus_slack_variables = 0;
    unsigned n_artificial_variables = 0;
    for(int i = 0; i < restrictionsTypes.size(); i++) {
        if(restrictionsTypes[i] == LESS_THAN_OR_EQUAL) n_surplus_slack_variables++;
        else if(restrictionsTypes[i] == GREATER_THAN_OR_EQUAL) {
            n_surplus_slack_variables++;
            n_artificial_variables++;
        }
        else if(restrictionsTypes[i] == EQUAL) n_artificial_variables++;
    }

    for(int i = 0; i < simplexTableau.columns(); i++) {
        zj.setElement(0, i, cb.dotProduct(simplexTableau.getColumn(i)));
    }
    cj_minus_zj = cj - zj;

    //simplexTableau.displayMatrix();
    //std::cout << std::endl;

    //b.displayMatrix();
    //std::cout << std::endl;
    //cj.displayMatrix();

    unsigned iterations = 0;
    while(!isSimplexDone(cj_minus_zj)) {

        //cj.displayMatrix();
        //std::cout << std::endl;

        pivots.setElement(0, 1, cj_minus_zj.maxValueIndex());
        Matrix ratios = zeros(restrictionsLHS.rows(), 1);

        ratios = b.pointDivision(simplexTableau.getColumn(pivots.getElement(0, 1)));
        //ratios.displayMatrix();
        //std::cout << std::endl;
        std::vector<double> simplexAux, bAux;
        for(int i = 0; i < ratios.rows(); i++) {
            simplexAux.push_back(simplexTableau.getElement(i, pivots.getElement(0, 1)));
            bAux.push_back(b.getElement(i, 0));
        }

        unsigned pivotRow = getPivotRow(simplexAux, bAux, ratios);
        if(pivotRow == -1) {
            optimalSolution = Matrix({INFINITY}, 1, 1);
            return Matrix({INFINITY}, 1, 1);
        }
        pivots.setElement(0, 0, pivotRow);

        unsigned oldBasis = pivots.getElement(0, 0);
        unsigned newBasis = pivots.getElement(0, 1);

        if(cb.getElement(oldBasis, 0) == M || cb.getElement(oldBasis, 0) == -1 * M) {
            unsigned artificial_index = basisIndices.getElement(oldBasis, 0);
            //std::cout << "art_index = " << artificial_index << std::endl;

            simplexTableau.removeColumn(artificial_index);
            cj.removeColumn(artificial_index);
            zj.removeColumn(artificial_index);
            cj_minus_zj.removeColumn(artificial_index);

            //if(newBasis > artificial_index) newBasis--;
            
            if(newBasis >= artificial_index) {
                //std::cout << "New basis out of bounds" << std::endl;
                newBasis -= 1;
            }

            for(int i = 0; i < basisIndices.rows(); i++) {
                if(basisIndices.getElement(i, 0) > artificial_index) basisIndices.setElement(i, 0, basisIndices.getElement(i, 0) - 1.0);
            }
        }

        basisIndices.setElement(oldBasis, 0, newBasis);
        cb.setElement(oldBasis, 0, cj.getElement(0, newBasis));

        //std::cout << "Pivot element is located at(" << oldBasis << ", " << newBasis - 1 << ") and is " << simplexTableau.getElement(oldBasis, newBasis - 1) << std::endl;
        //cb.displayMatrix();
        //std::cout << std::endl;
        //simplexTableau.displayMatrix();
        //std::cout << std::endl;

        Matrix newRow = Matrix(simplexTableau.getRow(oldBasis), 1, simplexTableau.columns()) * (1 / simplexTableau.getElement(oldBasis, newBasis));
        b.setElement(oldBasis, 0, b.getElement(oldBasis, 0) / simplexTableau.getElement(oldBasis, newBasis));

        //std::cout << "B after division: " << std::endl;
        //b.displayMatrix();
        //std::cout << std::endl << "--------------------------------------------------" << std::endl;

        //std::cout << "b factor = " << simplexTableau.getElement(oldBasis, newBasis) << std::endl;

        //std::cout << "B after division: " << std::endl;
        //b.displayMatrix();
        //std::cout << std::endl;

        //simplexTableau.displayMatrix();
        //std::cout << std::endl;
        //std::cout << b.getElement(oldBasis, 0) << std::endl;

        simplexTableau = simplexTableau.setRow(oldBasis, newRow);

        for(int i = 0; i < simplexTableau.rows(); i++) {
            if(i == oldBasis) continue;
            else {
                double factor = simplexTableau.getElement(i, newBasis);
                simplexTableau.rowOperation(oldBasis, i, -1 * factor);
                b.rowOperation(oldBasis, i, -1 * factor);
            }
        }

        simplexTableau.displayMatrix();
        std::cout << std::endl;
        b.displayMatrix();
        std::cout << std::endl;
        zj.displayMatrix();
        std::cout << std::endl;
        std::cout << "---------------------------------------------------------------" << std::endl;

        for(int i = 0; i < simplexTableau.columns(); i++) {
            zj.setElement(0, i, cb.dotProduct(simplexTableau.getColumn(i)));
        }
        cj_minus_zj = cj - zj;

        //cb.displayMatrix();
        //std::cout << std::endl;

        iterations++;

    }

    basisIndices.displayMatrix();
    std::cout << std::endl;

    simplexTableau.displayMatrix();
    std::cout << std::endl;
    zj.displayMatrix();
    std::cout << std::endl;
    cj_minus_zj.displayMatrix();
    std::cout << std::endl;
    //b.displayMatrix();
    //std::cout << std::endl;

    //cj.displayMatrix();
    //std::cout << std::endl;

    //basisIndices.displayMatrix();
    //std::cout << std::endl;

    Matrix solution = zeros(1, objectiveFunction.columns());

    //std::vector<double> solutionStuff(objectiveFunction.columns(), 0.0);


    for(int k = 0; k < basisIndices.rows(); k++) {
        if(basisIndices.getElement(k, 0) < objectiveFunction.columns()) solution.setElement(0, basisIndices.getElement(k, 0), b.getElement(k, 0));
        else if(basisIndices.getElement(k, 0) >= objectiveFunction.columns() + n_surplus_slack_variables && b.getElement(k, 0) > 0) {
            solution = Matrix({0}, 1, 1);
            break;
        }
    }

    //std::cout << "Hello there" << std::endl;
    optimalSolution = solution;
    return solution;
}

void LpProblem::displayProblem() {

    if(optimalSolution.rows() == 1 && optimalSolution.columns() == 1 && optimalSolution.getElement(0, 0) == 0) std::cout << "The problem is infeasible" << std::endl;
    else if(optimalSolution.rows() == 1 && optimalSolution.columns() == 1 && optimalSolution.getElement(0, 0) == INFINITY) std::cout << "The problem is unbounded" << std::endl;
    else {
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
        std::cout << "), and Z = " << std::setprecision(3) << std::fixed << optimalSolution.dotProduct(objectiveFunction) << std::endl;
    }
}

void LpProblem::addRestriction(Matrix lhs, restrictionType newType, double value) {
    restrictionsLHS.stackVertical(lhs);
    restrictionsRHS.stackVertical(Matrix({value}, 1, 1));
    restrictionsTypes.push_back(newType);
}

bool LpProblem::isProblemFeasible() {
    if(optimalSolution.rows() == 1 && optimalSolution.columns() == 1 && optimalSolution.getElement(0, 0) == 0) return false;
    else return true;
}

bool LpProblem::isProblemBounded() {
    if(optimalSolution.rows() == 1 && optimalSolution.columns() == 1 && optimalSolution.getElement(0, 0) == INFINITY) return false;
    else return true;
}

Matrix LpProblem::getOptimalSolution() {
    return optimalSolution;
}

ProblemType LpProblem::getType() {
    return type;
}

Matrix LpProblem::getObjectiveFunction() {
    return objectiveFunction;
}
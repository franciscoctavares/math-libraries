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

std::vector<std::pair<size_t, size_t>> LpProblem::getMaxWidth(std::vector<std::string> headers, Matrix matrix, std::string mode) {
    if(mode == "single column") {
        std::vector<std::pair<size_t, size_t>> maxWidth;
        size_t maxHeaderWidth = 0;
        size_t maxValueWidth = 0;

        for (size_t i = 0; i < headers.size(); ++i) {
            maxHeaderWidth = std::max(maxHeaderWidth, headers[i].length());
            maxValueWidth = std::max(maxValueWidth, std::to_string(matrix.getElement(i, 0)).length());
        }

        maxWidth.push_back(std::make_pair(maxHeaderWidth, maxValueWidth));
        return maxWidth;
    }
}

std::vector<std::string> LpProblem::basisHeaders(Matrix cj, Matrix basisIndices) {
    Matrix extraCj = cj.subMatrix(0, 0, objectiveFunction.columns(), cj.columns());
    Matrix basisIndexes = getBasisIndices(extraCj);
    basisIndexes = basisIndices;
    std::vector<std::vector<int>> constraintsIndexes = getRestrictionsIndexes(extraCj);

    std::vector<std::string> basisHeaders;
    for(int i = 0; i < basisIndexes.rows(); i++) {
        unsigned currentIndex = basisIndexes.getElement(i, 0);
        if(currentIndex >= objectiveFunction.columns()) {
            for(int k = 0; k < constraintsIndexes.size(); k++) {
                if(constraintsIndexes[k][0] == currentIndex) {
                    basisHeaders.push_back("s" + std::to_string(k + 1));
                    constraintsIndexes[k][0] = -2;
                    break;
                }
                else if(constraintsIndexes[k][1] == currentIndex) {
                    basisHeaders.push_back("a" + std::to_string(k + 1));
                    constraintsIndexes[k][1] = -2;
                    break;
                }
            }
        }
        else {
            //std::cout << "Its a variable x" << std::endl;
            basisHeaders.push_back("x" + std::to_string(currentIndex + 1));
        }
    }

    return basisHeaders;

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
    double total_z = cb.dotProduct(b);
    std::string total_z_str = std::to_string(total_z);
    total_z_str.erase(total_z_str.find_last_not_of('0') + 1, std::string::npos);
    if (total_z_str.back() == '.') total_z_str.pop_back(); // remove trailing dot if needed
    maxWidth = std::max(maxWidth, total_z_str.length());
    //maxWidth = std::max(maxWidth, std::to_string(cb.dotProduct(b)).length());

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

    std::ostringstream basis_oss;

    //std::cout << std::endl;



    std::vector<std::string> headers = {"x1", "x2", "x3", "s1", "s2", "s3", "s4", "s5", "s6"};
    std::vector<double> values = {30, 50, 40, 0, 0, 0, 0, 0, 0};

    std::vector<std::pair<size_t, size_t>> basisWidth = getMaxWidth(headers, cb, "single column");
    size_t maxHeaderWidth = basisWidth[0].first;
    size_t maxCoeffsWidth = basisWidth[0].second;

    maxHeaderWidth += 1;  // optional padding
    maxCoeffsWidth += 1;

    std::vector<std::string> basislines;

    size_t line_length_basis = maxHeaderWidth + maxCoeffsWidth + 2 + 2 + 2;
    line_length_basis--;
    
    // Spacing before cj row variables
    for(int i = 0; i < line_length_basis; i++) std::cout << " ";
    std::cout << "|";
    
    // variables row
    for(int i = 0; i < objectiveFunction.columns(); i++) {
        std::cout << std::setw(maxWidth) << "x" << i + 1 << " |";
    }
    for(int j = 0; j < restrictionsIndexes.size(); j++) {
        if(restrictionsIndexes[j][0] > 0) {
            std::cout << std::setw(maxWidth) << "s" << j + 1 << " |";
        }
    }
    for(int k = 0; k < restrictionsIndexes.size(); k++) {
        if(restrictionsIndexes[k][1] > 0) {

            std::cout << std::setw(maxWidth) << "a" << k + 1 << " |";
        }
    }

    std::cout << std::setw(maxWidth) << " " << "  |";
    
    // Elements of the cj row
    std::cout << std::endl;

    std::cout << "|" << std::setw(maxHeaderWidth) << "xB" << " |" << std::setw(maxCoeffsWidth + 1) << "cB" << " |";
    size_t line_length_basis_aux = line_length_basis - (5 + maxHeaderWidth);
    //for(int i = 0; i < line_length_basis_aux; i++) std::cout << " ";
    for(int i = 0; i < cj.columns(); i++) {
        std::cout << std::setw(maxWidth + 1) << cj.getElement(0, i) << " |";
    }
    std::cout << std::setw(maxWidth + 1) << "b" << " |" << std::endl;

    size_t total_line_length = line_length_basis + ((maxWidth + 3) * cj.columns()) + maxWidth + 1 + 2 + 1;
    for(size_t k = 0; k < total_line_length; k++) std::cout << "-";
    std::cout << std::endl;

    // Display basis variables and its objective function coefficients
    std::vector<std::string> headers_basis = basisHeaders(cj, basisIndexes);
    std::vector<double> basisCoeffs = cb.getElements();
    for(size_t i = 0; i < headers_basis.size(); i++) {
        std::cout << "|" << std::setw(maxHeaderWidth) << headers_basis[i] << " |"
                         << std::setw(maxCoeffsWidth + 1) << basisCoeffs[i] << " |";
        for(size_t j = 0; j < tableau.columns(); j++) {
            std::cout << std::setw(maxWidth + 1) << tableau.getElement(i, j) << " |";
        }

        std::cout << std::setw(maxWidth + 1) << b.getElement(i, 0) << " |";
        std::cout << std::endl;
    }

    for(size_t k = 0; k < total_line_length; k++) std::cout << "-";
    std::cout << std::endl;

    std::cout << std::setw(maxHeaderWidth + maxCoeffsWidth + 1 + 3) << "zj" << " |";

    for(int j = 0; j < zj.columns(); j++) {
        std::cout << std::setw(maxWidth + 1) << zj.getElement(0, j) << " |";
    }

    std::cout << std::setw(maxWidth + 1) << cb.dotProduct(b) << " |" << std::endl;

    std::cout << std::setw(maxHeaderWidth + maxCoeffsWidth + 1 + 3) << "cj - zj" << " |";

    for(int j = 0; j < zj.columns(); j++) {
        std::cout << std::setw(maxWidth + 1) << cj_minus_zj.getElement(0, j) << " |";
    }

    std::cout << std::endl << std::endl;
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
            else if(type == MAX) aux.push_back(-1 * M);
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
    //std::cout << "minValue = " << ratios.getElement(minIndex, 0) << " and minIndex = " << minIndex << std::endl << std::endl;
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

    //displaySimplexTableau(simplexTableau, cb, basisIndices, cj, b, zj, cj_minus_zj);

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

        simplexTableau = simplexTableau.setRow(oldBasis, newRow);

        for(int i = 0; i < simplexTableau.rows(); i++) {
            if(i == oldBasis) continue;
            else {
                double factor = simplexTableau.getElement(i, newBasis);
                simplexTableau.rowOperation(oldBasis, i, -1 * factor);
                b.rowOperation(oldBasis, i, -1 * factor);
            }
        }

        //simplexTableau.displayMatrix();
        //std::cout << std::endl;
        //b.displayMatrix();
        //std::cout << std::endl;
        //zj.displayMatrix();
        //std::cout << std::endl;
        //std::cout << "---------------------------------------------------------------" << std::endl;

        for(int i = 0; i < simplexTableau.columns(); i++) {
            zj.setElement(0, i, cb.dotProduct(simplexTableau.getColumn(i)));
        }
        cj_minus_zj = cj - zj;

        //displaySimplexTableau(simplexTableau, cb, basisIndices, cj, b, zj, cj_minus_zj);

        //std::cout << std::endl;
        //basisIndices.displayMatrix();
        //std::cout << std::endl;
        //cb.displayMatrix();
        //std::cout << std::endl;

        iterations++;

    }

    //displaySimplexTableau(simplexTableau, cb, basisIndices, cj, b, zj, cj_minus_zj);

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

void LpProblem::setOptimalSolution(Matrix newOptimalSolution) {
    optimalSolution = newOptimalSolution;
}

ProblemType LpProblem::getType() {
    return type;
}

Matrix LpProblem::getObjectiveFunction() {
    return objectiveFunction;
}

std::tuple<LpProblem, std::vector<std::pair<unsigned, unsigned>>, double> LpProblem::simplifyProblem(unsigned variable, double value) {
    double c = objectiveFunction.getElement(0, variable) * value;

    std::vector<std::pair<unsigned, unsigned>> oldVariablesToNewVariables;
    Matrix newObj = objectiveFunction;
    Matrix newLhs = restrictionsLHS;
    Matrix newRhs = restrictionsRHS;
    std::vector<restrictionType> newTypes = restrictionsTypes;

    newObj.removeColumn(variable);
    for(int i = 0; i < restrictionsLHS.rows(); i++) {
        double currentRhsValue = newRhs.getElement(i, 0);
        newRhs.setElement(i, 0, currentRhsValue - newLhs.getElement(i, variable) * value);
    }
    newLhs.removeColumn(variable);

    // verify that every rhs value is non negative. if not, multiply the entire constraint by -1
    for(int i = 0; i < newRhs.rows(); i++) {
        if(newRhs.getElement(i, 0) < 0) {
            newRhs.setElement(i, 0, newRhs.getElement(i, 0) * -1);
            if(newTypes[i] == LESS_THAN_OR_EQUAL) newTypes[i] = GREATER_THAN_OR_EQUAL;
            else if(newTypes[i] == GREATER_THAN_OR_EQUAL) newTypes[i] = LESS_THAN_OR_EQUAL;
            else newTypes[i] = EQUAL;
            for(int j = 0; j < newLhs.columns(); j++) {
                newLhs.setElement(i, j, newLhs.getElement(i, j) * -1);
            }
        }
    }

    // now, make the correspondences betwenn indices of the old variables and the new Variables
    if(variable == 0) {
        for(int j = 1; j < objectiveFunction.columns(); j++) {
            oldVariablesToNewVariables.push_back(std::make_pair(j, j - 1));
        }
    }
    else if(variable == objectiveFunction.columns() - 1) {
        for(int j = 0; j < objectiveFunction.columns() - 1; j++) {
            oldVariablesToNewVariables.push_back(std::make_pair(j, j));
        }
    }
    else {
        for(int j = 0; j < variable; j++) {
            oldVariablesToNewVariables.push_back(std::make_pair(j, j));
        }
        for(int j = variable + 1; j < objectiveFunction.columns(); j++) {
            oldVariablesToNewVariables.push_back(std::make_pair(j, j - 1));
        }     
    }

    LpProblem newProblem(type, newObj, newLhs, newRhs, newTypes);
    std::tuple<LpProblem, std::vector<std::pair<unsigned, unsigned>>, double> newProblemAndMoreStuff(newProblem, oldVariablesToNewVariables, c);

    return newProblemAndMoreStuff;
}

Matrix LpProblem::getSimplifiedProblemSolution(Matrix simplifiedSolution, std::vector<std::pair<unsigned, unsigned>> oldVarstoNewVars, unsigned variable, double value) {
    Matrix actualSolution = zeros(1, simplifiedSolution.columns() + 1);
    actualSolution.setElement(0, variable, value);

    for(int j = 0; j < oldVarstoNewVars.size(); j++) {
        double currentValue = simplifiedSolution.getElement(0, oldVarstoNewVars[j].second);
        actualSolution.setElement(0, oldVarstoNewVars[j].first, currentValue);
    }

    return actualSolution;
}

Matrix LpProblem::getConstraintsLHS() {
    return restrictionsLHS;
}

std::vector<restrictionType> LpProblem::getConstraintsTypes() {
    return restrictionsTypes;
}

Matrix LpProblem::getConstraintsRHS() {
    return restrictionsRHS;
}
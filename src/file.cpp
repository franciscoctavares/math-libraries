#include "../include/file.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <regex>
#include <string>
#include <algorithm>
#include <iomanip>
#include <cmath>

ModelFileReader::ModelFileReader() {
    
}

std::pair<std::string, std::vector<double>> ModelFileReader::readObjectiveFunction(std::string fileName) {
    std::ifstream file(fileName);
    if (!file) {
        std::cerr << "Error: Cannot open file.\n";
        //return 1;
    }

    std::string line;
    std::vector<double> objectiveCoefficients;
    std::vector<std::vector<double>> constraintsCoefficients;
    std::vector<double> rhsValues;
    std::vector<std::string> inequalitySigns;

    // --- Read objective function ---
    if (!std::getline(file, line)) {
        std::cerr << "Error: Missing objective function.\n";
        //return 1;
    }

    std::smatch match;
    std::regex headerRegex(R"((max|min):\s*(.*))");
    if (!std::regex_match(line, match, headerRegex)) {
        std::cerr << "Error: Invalid objective function format.\n";
        //return 1;
    }

    //std::cout << match[1] << std::endl;
    std::string problemType = match[1];

    std::string expr = match[2]; // part after 'max:' or 'min:'
    std::regex termRegex(R"(([+-]?\s*\d*\.?\d*)x(\d+))");

    int maxIndex = 0;
    std::vector<std::pair<int, double>> objTerms;

    auto termsBegin = std::sregex_iterator(expr.begin(), expr.end(), termRegex);
    auto termsEnd = std::sregex_iterator();

    for (auto it = termsBegin; it != termsEnd; ++it) {
        std::string coeffStr = (*it)[1].str();
        std::string indexStr = (*it)[2].str();

        coeffStr.erase(remove_if(coeffStr.begin(), coeffStr.end(), ::isspace), coeffStr.end());
        if (coeffStr.empty() || coeffStr == "+") coeffStr = "1";
        else if (coeffStr == "-") coeffStr = "-1";

        int varIndex = std::stoi(indexStr) - 1;
        double coeff = std::stod(coeffStr);
        objTerms.emplace_back(varIndex, coeff);
        maxIndex = std::max(maxIndex, varIndex);
    }

    objectiveCoefficients.resize(maxIndex + 1, 0.0);
    for (const auto& [index, value] : objTerms) {
        objectiveCoefficients[index] = value;
    }

    std::pair<std::string, std::vector<double>> result(problemType, objectiveCoefficients);
    return result;
}

std::vector<std::tuple<std::vector<double>, std::string, double>> ModelFileReader::readConstraints(std::string fileName) {
    std::ifstream file(fileName);
    if (!file) {
        std::cerr << "Error: Cannot open file.\n";
        //return 1;
    }

    std::string line;
    std::vector<double> objectiveCoefficients;
    std::vector<std::vector<double>> constraintsCoefficients;
    std::vector<double> rhsValues;
    std::vector<std::string> inequalitySigns;

    std::vector<std::tuple<std::vector<double>, std::string, double>> allConstraints;

    // --- Read objective function ---
    if (!std::getline(file, line)) {
        std::cerr << "Error: Missing objective function.\n";
        //return 1;
    }
    std::regex termRegex(R"(([+-]?\s*\d*\.?\d*)x(\d+))");
    /*
    std::smatch match;
    std::regex headerRegex(R"((max|min):\s*(.*))");
    if (!std::regex_match(line, match, headerRegex)) {
        std::cerr << "Error: Invalid objective function format.\n";
        //return 1;
    }

    //std::cout << match[1] << std::endl;

    std::string expr = match[2]; // part after 'max:' or 'min:'
    std::regex termRegex(R"(([+-]?\s*\d*\.?\d*)x(\d+))");

    int maxIndex = 0;
    std::vector<std::pair<int, double>> objTerms;

    auto termsBegin = std::sregex_iterator(expr.begin(), expr.end(), termRegex);
    auto termsEnd = std::sregex_iterator();

    for (auto it = termsBegin; it != termsEnd; ++it) {
        std::string coeffStr = (*it)[1].str();
        std::string indexStr = (*it)[2].str();

        coeffStr.erase(remove_if(coeffStr.begin(), coeffStr.end(), ::isspace), coeffStr.end());
        if (coeffStr.empty() || coeffStr == "+") coeffStr = "1";
        else if (coeffStr == "-") coeffStr = "-1";

        int varIndex = std::stoi(indexStr) - 1;
        double coeff = std::stod(coeffStr);
        objTerms.emplace_back(varIndex, coeff);
        maxIndex = std::max(maxIndex, varIndex);
    }

    objectiveCoefficients.resize(maxIndex + 1, 0.0);
    for (const auto& [index, value] : objTerms) {
        objectiveCoefficients[index] = value;
    }
    */

    // --- Skip blank line ---
    while (std::getline(file, line)) {
        if (line.empty()) break;
    }

    // --- Read constraints ---
    std::regex constraintRegex(R"(.*?(<=|>=|=)\s*(-?\d+\.?\d*))");

    while (std::getline(file, line)) {
        if (line.empty()) continue;

        std::smatch rhsMatch;
        if (!std::regex_search(line, rhsMatch, constraintRegex)) {
            std::cerr << "Warning: Skipping invalid constraint line: " << line << "\n";
            continue;
        }

        std::string inequality = rhsMatch[1];
        double rhs = std::stod(rhsMatch[2]);
        inequalitySigns.push_back(inequality);
        rhsValues.push_back(rhs);

        std::string lhsPart = line.substr(0, rhsMatch.position(1)); // LHS up to inequality
        std::vector<std::pair<int, double>> constraintTerms;
        int maxVarInConstraint = 0;

        auto lhsTermsBegin = std::sregex_iterator(lhsPart.begin(), lhsPart.end(), termRegex);
        auto lhsTermsEnd = std::sregex_iterator();

        for (auto it = lhsTermsBegin; it != lhsTermsEnd; ++it) {
            std::string coeffStr = (*it)[1].str();
            std::string indexStr = (*it)[2].str();

            coeffStr.erase(remove_if(coeffStr.begin(), coeffStr.end(), ::isspace), coeffStr.end());
            if (coeffStr.empty() || coeffStr == "+") coeffStr = "1";
            else if (coeffStr == "-") coeffStr = "-1";

            int varIndex = std::stoi(indexStr) - 1;
            double coeff = std::stod(coeffStr);
            constraintTerms.emplace_back(varIndex, coeff);
            maxVarInConstraint = std::max(maxVarInConstraint, varIndex);
        }

        std::vector<double> coeffs(std::max((int)objectiveCoefficients.size(), maxVarInConstraint + 1), 0.0);
        for (const auto& [index, value] : constraintTerms) {
            coeffs[index] = value;
        }

        constraintsCoefficients.push_back(coeffs);

        std::tuple<std::vector<double>, std::string, double> constraint(coeffs, inequality, rhs);
        allConstraints.push_back(constraint);
    }

    return allConstraints;
}

std::pair<std::pair<std::string, std::vector<double>>, std::vector<std::tuple<std::vector<double>, std::string, double>>> ModelFileReader::readModel(std::string fileName) {
    std::pair<std::string, std::vector<double>> objectiveFunction = readObjectiveFunction(fileName);
    std::vector<std::tuple<std::vector<double>, std::string, double>> constraints = readConstraints(fileName);
    return std::make_pair(objectiveFunction, constraints);
}

void ModelFileReader::displayModel(std::pair<std::pair<std::string, std::vector<double>>, std::vector<std::tuple<std::vector<double>, std::string, double>>> modelParts) {
    std::pair<std::string, std::vector<double>> objectiveFunction = modelParts.first;
    std::vector<std::tuple<std::vector<double>, std::string, double>> constraints = modelParts.second;

    // Display objective function
    std::cout << objectiveFunction.first << " Z = ";
    for(int i = 0; i < objectiveFunction.second.size(); i++) {
        double coeff = objectiveFunction.second[i];
        if(coeff < 0) std::cout << "- ";
        else { if(i > 0) std::cout << "+ "; }

        if(fabs(coeff) != 1) std::cout << fabs(coeff);
        std::cout << "x" << i + 1;
        if(i < objectiveFunction.second.size() - 1) std::cout << " ";
    }
    std::cout << std::endl << std::endl;

    std::cout << "Subject to:" << std::endl << std::endl;

    // Display constraints
    for(std::tuple<std::vector<double>, std::string, double> constraint: constraints) {
        std::vector<double> lhs = std::get<0>(constraint);
        std::string restrictionType = std::get<1>(constraint);
        double rhs = std::get<2>(constraint);

        bool hasWrittenCoeffs = false;
        for(int i = 0; i < lhs.size(); i++) {
            if(lhs[i] == 0) continue;
            double coeff = lhs[i];
            if(coeff < 0) std::cout << "- ";
            else { if(hasWrittenCoeffs) std::cout << "+ "; }

            if(fabs(coeff) != 1) std::cout << fabs(coeff);
            hasWrittenCoeffs = true;
            std::cout << "x" << i + 1;
            if(i < lhs.size() - 1) std::cout << " ";
        }

        std::cout << " " << restrictionType << " " << rhs << std::endl;
    }

}
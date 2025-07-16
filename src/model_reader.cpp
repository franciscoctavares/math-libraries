#include "../include/model_reader.h"

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

std::vector<std::tuple<std::vector<double>, std::string, double>> ModelFileReader::readConstraints(std::string fileName, unsigned varsNumber) {
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
        int globalMaxVarIndex = 0;

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
            globalMaxVarIndex = std::max(globalMaxVarIndex, varIndex);
        }

        std::vector<double> coeffs(varsNumber, 0.0);
        for (const auto& [index, value] : constraintTerms) {
            coeffs[index] = value;
        }

        constraintsCoefficients.push_back(coeffs);

        //constraintsCoefficients[0].resize();

        std::tuple<std::vector<double>, std::string, double> constraint(coeffs, inequality, rhs);
        allConstraints.push_back(constraint);
    }

    return allConstraints;
}

LpProblem ModelFileReader::readModel(std::string fileName) {
    std::pair<std::string, std::vector<double>> objectiveFunctionAux = readObjectiveFunction(fileName);
    std::vector<std::tuple<std::vector<double>, std::string, double>> constraintsAux = readConstraints(fileName, objectiveFunctionAux.second.size());

    ProblemType type;
    if(objectiveFunctionAux.first == "max") type = MAX;
    else if(objectiveFunctionAux.first == "min") type = MIN;

    std::vector<double> objectiveFunction = objectiveFunctionAux.second;
    std::vector<Constraint> constraints;
    for(std::tuple<std::vector<double>, std::string, double> consts : constraintsAux) {
        constraints.push_back(Constraint(std::get<0>(consts), std::get<1>(consts),std::get<2>(consts)));
    }

    return LpProblem(type, objectiveFunction, constraints);
}
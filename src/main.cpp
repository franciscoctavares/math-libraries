#include <iostream>
#include <vector>
#include "include/matrix.h"
#include "include/constraint.h"
#include "include/lp.h"
#include "include/file.h"

using namespace std;

int main() {
    /*
    ProblemType probType = MAX;
    vector<double> objectiveFunction = {30, 50, 40};
    Constraint constraint1({1, 2, 2}, "<=", 100);
    Constraint constraint2({2, 1, 2}, "<=", 80);
    Constraint constraint3({2, 2, 1}, "<=", 60);
    vector<Constraint> constraints = {constraint1, constraint2, constraint3};

    LpProblem problem(probType, objectiveFunction, constraints);

    Matrix solution = problem.solveSimplex();
    problem.displayProblem(solution);
    */

    ModelFileReader reader;
    std::pair<std::pair<std::string, std::vector<double>>, std::vector<std::tuple<std::vector<double>, std::string, double>>> model = reader.readModel("model.lp");
    
    ProblemType type;
    if(model.first.first == "max") type = MAX;
    else if(model.first.first == "min") type = MIN;

    std::vector<double> objectiveFunction = model.first.second;
    std::vector<Constraint> constraints;
    for(std::tuple<std::vector<double>, std::string, double> consts : model.second) {
        constraints.push_back(Constraint(std::get<0>(consts), std::get<1>(consts),std::get<2>(consts)));
    }

    LpProblem problem2(type, objectiveFunction, constraints);
    problem2.solveSimplex();
    problem2.displayProblem(problem2.getOptimalSolution());

    return 0;
}
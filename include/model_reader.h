#ifndef FILE_H
#define FILE_H

#include <vector>
#include <string>
#include <tuple>
#include <utility>
#include "lp.h"


class ModelFileReader {
    private:
        std::pair<std::string, std::vector<double>> readObjectiveFunction(std::string);
        std::vector<std::tuple<std::vector<double>, std::string, double>> readConstraints(std::string, unsigned varsNumber);
    public:
        ModelFileReader();
        
        /**
         * @brief Given the path of the model file, reads a linear programming problem and returns an LpProblem object with the model
         * 
         * @param fileName the path to the model file
         * @return LpProblem - the model read, encapsulated in the LpProblem class
         */
        LpProblem readModel(std::string fileName);
};

#endif
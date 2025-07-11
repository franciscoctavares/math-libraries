#ifndef FILE_H
#define FILE_H

#include <vector>
#include <string>
#include <tuple>
#include <utility>


class ModelFileReader {
    public:
        ModelFileReader();
        std::pair<std::string, std::vector<double>> readObjectiveFunction(std::string);
        std::vector<std::tuple<std::vector<double>, std::string, double>> readConstraints(std::string);
        std::pair<std::pair<std::string, std::vector<double>>, std::vector<std::tuple<std::vector<double>, std::string, double>>> readModel(std::string);
        void displayModel(std::pair<std::pair<std::string, std::vector<double>>, std::vector<std::tuple<std::vector<double>, std::string, double>>>);
};

#endif
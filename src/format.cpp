#include "../include/format.h"

std::vector<std::pair<size_t, size_t>> getMaxWidth(std::vector<std::string> headers, std::vector<std::vector<double>> matrix, std::string mode) {
    if(mode == "column") {
        std::vector<std::pair<size_t, size_t>> maxWidth;
        size_t maxHeaderWidth = 0;
        size_t maxValueWidth = 0;

        for (size_t i = 0; i < headers.size(); ++i) {
            maxHeaderWidth = std::max(maxHeaderWidth, headers[i].length());
            maxValueWidth = std::max(maxValueWidth, std::to_string(matrix[i][0]).length());
        }

        maxWidth.push_back(std::make_pair(maxHeaderWidth, maxValueWidth));
        return maxWidth;
    }
    else if(mode == "matrix") {
        // in matrix mode, computes the max width per row and then returns the max width of all the rows
        
    }
}
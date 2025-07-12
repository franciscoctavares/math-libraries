#ifndef FORMAT_H
#define FORMAT_H

#include <cstddef>
#include <vector>
#include <string>

std::vector<std::pair<size_t, size_t>> getMaxWidth(std::vector<std::string>, std::vector<std::vector<double>>, std::string mode = "column");

#endif

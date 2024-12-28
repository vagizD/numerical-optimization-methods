#pragma once
#include "datasets.h"

namespace ADAAI {
    double computeSortinoRatio(const std::vector<double>& params, const Dataset& dataset);
    double computePnL(const std::vector<double>& params, const Dataset& data);
}

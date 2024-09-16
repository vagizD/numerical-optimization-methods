#pragma once
#include <tuple>
#include "model_params.h"

std::tuple<std::size_t, std::size_t, double, double> adjust_grid_boundary(
    const ModelParams &mp,
    std::size_t n,
    std::size_t m,
    double T,
    double M,
    double S_0,
    double V_0,
    double K
);

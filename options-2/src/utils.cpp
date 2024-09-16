#include "utils.h"
#include <cmath>
#include <iostream>
#include "model_params.h"

// computes S_max, V_max estimates and adjusts them such that S_0 and V_0 lie on a grid
// S_max and V_max are both increased such that a = S_0 * n / S_max becomes the next
// closest integer, i.e. i* = [a]
std::tuple<size_t, size_t, double, double> adjust_grid_boundary(
    const ModelParams &mp,
    const size_t n,
    const size_t m,
    const double T,
    const double M,
    const double S_0,
    const double V_0,
    const double K
) {
    // Rule of thumb
    double S_max = std::max(S_0, K) * std::exp(M * mp.sigma * std::sqrt(T));
    double V_max = std::max(V_0, mp.theta) * (1 + M * mp.epsilon * std::sqrt(T));

    double i_double = S_0 * static_cast<int>(n) / S_max;
    double j_double = V_0 * static_cast<int>(m) / V_max;

    double S_delta = S_max * (i_double / std::round(i_double) - 1.0);
    double V_delta = V_max * (j_double / std::round(j_double) - 1.0);

    auto i = static_cast<size_t>(std::round(i_double));
    auto j = static_cast<size_t>(std::round(j_double));

    return std::make_tuple(i, j, S_max + S_delta, V_max + V_delta);
}

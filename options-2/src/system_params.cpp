#include "system_params.h"
#include <cmath>
#include <iostream>

std::tuple<size_t, size_t, double, double> SystemParams::adjust_grid_boundary(
    const ModelParams &a_mp,
    const size_t a_n,
    const size_t a_m,
    const double a_T,
    const double a_M,
    const double a_S_0,
    const double a_V_0,
    const double a_K
) {
    // Rule of thumb
    const double S_max =
        std::max(a_S_0, a_K) * std::exp(a_M * a_mp.m_sigma * std::sqrt(a_T));
    const double V_max =
        std::max(a_V_0, a_mp.m_theta) * (1 + a_M * a_mp.m_epsilon * std::sqrt(a_T));

    const double j_double = a_S_0 * static_cast<int>(a_n) / S_max;
    const double i_double = a_V_0 * static_cast<int>(a_m) / V_max;

    const double S_delta = S_max * (j_double / std::round(j_double) - 1.0);
    const double V_delta = V_max * (i_double / std::round(i_double) - 1.0);

    const auto j = static_cast<size_t>(std::round(j_double));
    const auto i = static_cast<size_t>(std::round(i_double));

    return std::make_tuple(j, i, S_max + S_delta, V_max + V_delta);
}

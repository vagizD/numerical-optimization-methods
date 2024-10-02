#pragma once

#pragma once

#include <cmath>
#include "gsl/gsl_cdf.h"
#include "gsl/gsl_randist.h"
#include "system_params.h"

double option_price(const SystemParams &a_sp) {
    double V = a_sp.m_mp.m_sigma;
    double int_risk = a_sp.m_mp.m_r;
    double dp = (std::log(a_sp.m_S_max / (double)a_sp.m_K) + int_risk + V / 2) /
                std::sqrt(V);  // d+
    double dm = (std::log(a_sp.m_S_max / (double)a_sp.m_K) + int_risk - V / 2) /
                std::sqrt(V);  // d-

    // s_max * Ф(d+) - K * e^{-integrate_risk} * Ф(d-)
    return a_sp.m_S_max * gsl_cdf_ugaussian_P(dp) -
           a_sp.m_K * std::exp(-int_risk) * gsl_cdf_ugaussian_P(dm);
}

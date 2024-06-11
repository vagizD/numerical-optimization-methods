#pragma once

#include "constants.h"
#include "gsl/gsl_cdf.h"
#include "gsl/gsl_randist.h"

double option_price(double s_max, double t_max) {
    double V = integrate_volatility(t_max);
    double int_risk = integrate_risk(t_max);
    double dp = (std::log(s_max / (double)K) + int_risk + V / 2) / std::sqrt(V);  // d+
    double dm = (std::log(s_max / (double)K) + int_risk - V / 2) / std::sqrt(V);  // d-

    // s_max * Ф(d+) - K * e^{-integrate_risk} * Ф(d-)
    return s_max * gsl_cdf_ugaussian_P(dp) -
           K * std::exp(-int_risk) * gsl_cdf_ugaussian_P(dm);
}

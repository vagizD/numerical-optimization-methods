#pragma once

#include <array>
#include <cassert>
#include <cinttypes>
#include <complex>

constexpr size_t N = 500;      // Number of discrete steps for S (price of share)
constexpr int32_t K = 100;     // GazProm share price
constexpr double STEP = 1e-3;  // Integration step

double risk(double tau) {
    assert(0 <= tau);
    assert(tau <= 1);
    if (tau <= 0.25) {
        return 0.16;
    }
    if (tau <= 0.50) {
        return 0.15;
    }
    if (tau <= 0.75) {
        return 0.14;
    }
    return 0.13;
}

double volatility(double tau) {
    assert(0 <= tau);
    assert(tau <= 1);
    if (tau <= 0.25) {
        return 0.23;
    }
    if (tau <= 0.50) {
        return 0.27;
    }
    if (tau <= 0.75) {
        return 0.25;
    }
    return 0.28;
}

double integrate_risk(double tau) {
    double res = 0;
    if (tau <= 0.25) {
        return res + tau * risk(0);
    }
    res += 0.25 * risk(0);
    if (tau <= 0.5) {
        return res + (tau - 0.25) * risk(0.30);
    }
    res += 0.25 * risk(0.30);
    if (tau <= 0.75) {
        return res + (tau - 0.50) * risk(0.60);
    }
    res += 0.25 * risk(0.6);
    return res + (tau - 0.75) * risk(1);
}

// S_max = K * exp(N * volatility_max * tau_max^(1/2))
const double S_max = K * std::exp(4 * 0.28);

constexpr std::array<double, 1> k2_consts = {1. / 4.};
constexpr std::array<double, 2> k3_consts = {3. / 32., 9. / 32.};
constexpr std::array<double, 3> k4_consts = {
    1932. / 2197., -7200. / 2197., 7296. / 2197.};
constexpr std::array<double, 4> k5_consts = {
    439. / 216., -8., 3680. / 513., -845. / 4104.};
constexpr std::array<double, 5> k6_consts = {
    -8. / 27., 2., -3544. / 2565., 1859. / 4104., -11. / 40.};
constexpr std::array<double, 6> delta_consts = {-1. / 360.,     0.,        128. / 4275.,
                                                2197. / 75240., -1. / 50., -2. / 55.};
constexpr std::array<double, 6> gamma_consts = {
    16. / 135, 0., 6656. / 12825., 28561. / 56430., -9. / 50., 2. / 55.};

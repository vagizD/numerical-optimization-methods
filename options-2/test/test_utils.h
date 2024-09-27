#pragma once
#include <cassert>
#include <iomanip>
#include <iostream>
#include "constants.h"
#include "utils.h"

inline void test_adjust_boundary() {
    std::cout << "TEST: ADJUST BOUNDARY" << '\n';
    constexpr double r = 0.25;
    constexpr double d = 0.20;
    constexpr double sigma = 0.16;
    constexpr double kappa = 0.10;
    constexpr double beta = 1.0;
    constexpr double epsilon = 0.13;
    constexpr double rho = 0.50;
    constexpr double theta = 0.30;

    constexpr size_t n = 200;
    constexpr size_t m = 200;
    constexpr size_t p = 100;
    constexpr double T = 1.0;
    constexpr double M = 1.0;
    constexpr double S_0 = 100.0;
    constexpr double V_0 = 1.0;
    constexpr double K = 90;

    auto sp = SystemParams(
        r, d, sigma, kappa, beta, epsilon, rho, theta, S_0, V_0, T, M, K, n, m, p
    );

    const double S_max_est =
        S_0 * static_cast<double>(n) / static_cast<double>(sp.m_S_0_j);
    const double V_max_est =
        V_0 * static_cast<double>(m) / static_cast<double>(sp.m_V_0_i);
    std::cout << std::setprecision(5) << std::fixed;
    std::cout << "i          : " << sp.m_V_0_i << '\n';
    std::cout << "j          : " << sp.m_S_0_j << '\n';
    std::cout << "S_max      : " << sp.m_S_max << '\n';
    std::cout << "S_0 * n / j: " << S_max_est << '\n';
    std::cout << "V_max      :   " << sp.m_V_max << '\n';
    std::cout << "V_0 * m / i:   " << V_max_est << '\n';

    assert(std::abs(S_max_est - sp.m_S_max) <= RTOL);
    assert(std::abs(V_max_est - sp.m_V_max) <= RTOL);

    std::cout << "PASSED\n\n";
}

inline void test_utils() {
    std::cout << "\n=================== START TEST UTILS ===================\n";
    test_adjust_boundary();
    std::cout << "==================== END TEST SLAE ====================\n\n";
}
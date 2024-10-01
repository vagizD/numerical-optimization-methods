#pragma once
#include <array>
#include <iostream>
#include <chrono>
#include <vector>
#include "container.h"
#include "stepper.h"

inline void test_stepper() {
    std::cout << "\n=================== START TEST STEPPER ===================\n";
    auto t1 = std::chrono::high_resolution_clock::now();

    constexpr double r = 0.19;
    constexpr double d = 0.01;
    constexpr double sigma = 0.20;
    constexpr double kappa = 0.10;
    constexpr double beta = 0.90;
    constexpr double epsilon = 0.30;
    constexpr double rho = -0.15;
    constexpr double theta = 1.0;

    constexpr size_t n = 20;
    constexpr size_t m = 20;
    constexpr size_t p = 365;
    constexpr double T = 1.0;
    constexpr double M = 1.0;
    constexpr double S_0 = 100.0;
    constexpr double V_0 = 1.0;
    constexpr double K = 105;

    const auto sp = SystemParams(
        r, d, sigma, kappa, beta, epsilon, rho, theta, S_0, V_0, T, M, K, n, m, p
    );

    std::vector<std::array<double, (n + 1) * (m + 1)>> results;

    {
        Matrix2D<(n + 1) * (m + 1)> A;
        std::array<double, (n + 1) * (m + 1)> res = {};
        option_price_stepper<Matrix2D<(n + 1) * (m + 1)>, m, n>(
            sp, A, res
        );

        results.push_back(res);
    }

    {
        Matrix1D<(n + 1) * (m + 1)> A;
        std::array<double, (n + 1) * (m + 1)> res = {};
        option_price_stepper<Matrix1D<(n + 1) * (m + 1)>,  m, n>(
            sp, A, res
        );

        results.push_back(res);
    }


    {
        MatrixGSL<(n + 1) * (m + 1)> A;
        std::array<double, (n + 1) * (m + 1)> res = {};
        option_price_stepper<
            MatrixGSL<(n + 1) * (m + 1)>, m, n>(sp, A, res);

        results.push_back(res);
    }

    std::cout << "    (V_i, S_j)    == Matrix2D == Matrix1D == MatrixGSL \n";

    for (int i = 0; i < m + 1; i++) {
        for (int j = 0; j < n + 1; j++) {
            std::cout << "(" << sp.m_V_max * i / m << ", " << sp.m_S_max * j / n
                      << ")   ";
            for (int l = 0; l < results.size(); l++) {
                std::cout << results[l][i * (n + 1) + j] << "  | ";
            }
            std::cout << "\n";
        }
    }
    auto t2 = std::chrono::high_resolution_clock::now();
    auto ms_int = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
    std::cout << static_cast<double>(ms_int.count()) / 1e3 << " seconds\n";
    std::cout << "=================== END TEST STEPPER ===================\n";
}

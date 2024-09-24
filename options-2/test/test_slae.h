#pragma once
#include <cassert>
#include <iomanip>
#include <iostream>
#include "slae.h"
#include "test_constants.h"

template <const size_t N>
void run_test(
    const std::array<double, N * N> &A,
    const std::array<double, N> &b,
    const std::array<double, N> &ans,
    const size_t test_id
) {
    std::array<double, N> x{};

    solve_slae(A, b, x);

    std::cout << "TEST " << test_id << '\n';
    for (size_t i = 0; i < N; ++i) {
        std::cout << "x" << i + 1 << " = " << x[i] << ", ans" << i + 1 << " = " << ans[i]
                  << '\n';
        assert(std::abs(x[i] - ans[i]) <= RTOL);
    }
    std::cout << "PASSED\n\n";
}

inline void test_slae() {
    std::cout << "\n=================== START TEST SLAE ===================\n";
    constexpr size_t N1 = 3;
    constexpr std::array<double, N1 *N1> A1 = {1, 0, 0, 0, 1, 0, 0, 0, 1};
    constexpr std::array<double, N1> b1 = {1, 2, 3};
    constexpr std::array<double, N1> ans1 = {1, 2, 3};
    run_test<N1>(A1, b1, ans1, 1);

    constexpr size_t N2 = 3;
    constexpr std::array<double, N2 *N2> A2 = {1, 1, 0, 0, 1, 1, 1, 0, 1};
    constexpr std::array<double, N2> b2 = {1, 2, 3};
    constexpr std::array<double, N1> ans2 = {1, 0, 2};
    run_test<N2>(A2, b2, ans2, 2);

    std::cout << "\n==================== END TEST SLAE ====================\n\n";
}

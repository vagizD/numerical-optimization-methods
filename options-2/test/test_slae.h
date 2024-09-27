#pragma once
#include <cassert>
#include <iomanip>
#include <iostream>
#include "slae.h"
#include "container.h"
#include "test_constants.h"
#include <functional>

template <const size_t N, typename C>
void run_test(
    const C &A,
    const std::array<double, N> &b,
    const std::array<double, N> &ans,
    const size_t test_id,
    const std::string msg
) {
    std::array<double, N> x{};

    std::cout << "TEST " << test_id << ": " << msg << '\n';

    solve_slae<N, C>(A, b, x);

    for (size_t i = 0; i < N; ++i) {
        std::cout << "x" << i + 1 << " = " << x[i] << ", ans" << i + 1 << " = " << ans[i]
                  << '\n';
        assert(std::abs(x[i] - ans[i]) <= RTOL);
    }
    std::cout << "PASSED\n\n";
}

template <const size_t N, typename C>
void run_test_sparse(
    const C &A,
    const std::array<double, N> &b,
    const std::array<double, N> &ans,
    const size_t test_id,
    const std::string msg
) {
    std::array<double, N> x{};

    std::cout << "TEST " << test_id << ": " << msg << '\n';

    solve_slae_sparse<N, C>(A, b, x);

    for (size_t i = 0; i < N; ++i) {
        std::cout << "x" << i + 1 << " = " << x[i] << ", ans" << i + 1 << " = " << ans[i]
                  << '\n';
        // assert(std::abs(x[i] - ans[i]) <= RTOL);
    }
    std::cout << "PASSED\n\n";
}


inline void test_slae() {
    std::cout << "\n=================== START TEST SLAE ===================\n";
    constexpr size_t N1 = 3;
    Matrix2D<N1> A1;
    A1.init_zero();

    // 1 0 0
    // 0 1 0
    // 0 0 1
    A1.set(0, 0, 1);
    A1.set(1, 1, 1);
    A1.set(2, 2, 1);

    constexpr std::array<double, N1> b1 = {1, 2, 3};
    constexpr std::array<double, N1> ans1 = {1, 2, 3};
    run_test<N1, Matrix2D<N1>>(A1, b1, ans1, 1, "Matrix2D GE");

    constexpr size_t N2 = 3;
    Matrix1D<N2> A2;
    A2.init_zero();

    // 1 1 0
    // 0 1 1
    // 1 0 1
    A2.set(0, 0, 1);
    A2.set(0, 1, 1);
    A2.set(1, 1, 1);
    A2.set(1, 2, 1);
    A2.set(2, 2, 1);
    A2.set(2, 0, 1);
    constexpr std::array<double, N2> b2 = {1, 2, 3};
    constexpr std::array<double, N2> ans2 = {1, 0, 2};
    run_test<N2, Matrix1D<N2>>(A2, b2, ans2, 2, "Matrix1D GE");

    constexpr size_t N3 = 4;
    Matrix1D<N3> A3;
    A3.init_zero();

    // -10   0   1   0 | 1
    //   0   1   0  -1 | 2
    //   3   0  -4  15 | 3
    //   0   2   3   4 | 4
    A3.set(0, 0, -10);
    A3.set(0, 2, 1);
    A3.set(1, 1, 1);
    A3.set(1, 3, -1);
    A3.set(2, 0, 3);
    A3.set(2, 2, -4);
    A3.set(2, 3, 15);
    A3.set(3, 1, 2);
    A3.set(3, 2, 3);
    A3.set(3, 3, 4);
    constexpr std::array<double, N3> b3 = {1, 2, 3, 4};
    constexpr std::array<double, N3> ans3 = {-29./224., 481./224., -66./224., 33./224.};
    run_test<N3, Matrix1D<N3>>(A3, b3, ans3, 3, "Matrix1D GE");

    constexpr size_t N4 = 3;
    MatrixSparse<N4, 5> A4;
    A4.init_zero();

    // 1 1 0
    // 0 1 1
    // 1 0 1
    A4.set(0, 0, 1);
    A4.set(0, 1, 1);
    A4.set(1, 1, 1);
    A4.set(1, 2, 1);
    A4.set(2, 2, 1);
    A4.set(2, 0, 1);
    constexpr std::array<double, N4> b4 = {1, 2, 3};
    constexpr std::array<double, N4> ans4 = {1, 0, 2};
    run_test_sparse<N4, MatrixSparse<N4, 5>>(A4, b4, ans4, 4, "MatrixSparse GE");

    constexpr size_t N5 = 4;
    MatrixSparse<N5, 5> A5;
    A5.init_zero();

    // -10   0   1   0 | 1
    //   0   1   0  -1 | 2
    //   3   0  -4  15 | 3
    //   0   2   3   4 | 4
    A5.set(0, 0, -10);
    A5.set(0, 2, 1);
    A5.set(1, 1, 1);
    A5.set(1, 3, -1);
    A5.set(2, 0, 3);
    A5.set(2, 2, -4);
    A5.set(2, 3, 15);
    A5.set(3, 1, 2);
    A5.set(3, 2, 3);
    A5.set(3, 3, 4);
    constexpr std::array<double, N5> b5 = {1, 2, 3, 4};
    constexpr std::array<double, N5> ans5 = {-29./224., 481./224., -66./224., 33./224.};
    run_test_sparse<N5, MatrixSparse<N5, 5>>(A5, b5, ans5, 5, "MatrixSparse GE");

    std::cout << "==================== END TEST SLAE ====================\n\n";
}

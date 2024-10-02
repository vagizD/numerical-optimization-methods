#pragma once
#include <cassert>
#include <iomanip>
#include <iostream>
#include <fstream>
#include "constants.h"
#include "container.h"
#include "slae.h"
#include "blas.h"

template <typename F_ge, typename F_gsl, const size_t N>
void save(
    const std::array<size_t, N> &Ns,
    const std::array<F_ge, N> &norms_ge,
    const std::array<F_gsl, N> &norms_gsl,
    const std::string &filename
) {
    std::cout << "SAVING...\n";
    std::ofstream file;
    file.open(filename);
    file << "n,norm_ge,norm_gsl\n";
    for (size_t i = 0; i < N; ++i) {
        file << Ns[i] << ',' << norms_ge[i] << ',' << norms_gsl[i] << '\n';
    }
    file.close();
}

template <typename F, const size_t N, typename C>
void run_test_with_ans(
    const C &A,
    const std::array<F, N> &b,
    const std::array<F, N> &ans,
    const size_t test_id,
    const std::string &msg
) {
    std::array<F, N> x{};

    std::cout << "TEST " << test_id << ": " << msg << '\n';

    solve_slae<F, N, C>(A, b, x);

    for (size_t i = 0; i < N; ++i) {
        std::cout << "x" << i + 1 << " = " << x[i] << ", ans" << i + 1 << " = " << ans[i]
                  << '\n';
        assert(std::abs(x[i] - ans[i]) <= RTOL);
    }
    std::cout << "PASSED\n\n";
}

template <typename F, const size_t N, typename C>
F run_test_no_ans(
    const C &A,
    const std::array<F, N> &b,
    const size_t test_id,
    const std::string &msg,
    const bool verbose
) {
    std::array<F, N> x{};
    std::array<F, N> Ax;
    std::array<F, N> diff{};

    if (verbose)
        std::cout << "TEST " << test_id << ": " << msg << '\n';

    solve_slae<F, N, C>(A, b, x);
    sqmatrix_on_vector<F, C, N>(A, x, Ax);
    vectors_diff<F, N>(Ax, b, diff);
    F res_norm = vector_norm<F, N>(diff);

    if (verbose)
        std::cout << "N = " << N <<  ", ||Ax - b|| = " << res_norm << '\n';
//    assert(res_norm <= RTOL);
    if (verbose)
        std::cout << "PASSED\n\n";
    return res_norm;
}

template<typename F, const size_t N, typename C>
F make_normal_test(int test_id) {
    std::cout << "NORMAL TEST: " << test_id << '\n';
    C A;
    std::array<F, N> b{};
    A.init_normal();
    generate_normal<F, N>(b, 0.0, 1.0);
    return run_test_no_ans(A, b, -1, "", false);
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
    run_test_with_ans<double, N1, Matrix2D<N1>>(A1, b1, ans1, 1, "Matrix2D GE");

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
    run_test_with_ans<double, N2, Matrix1D<N2>>(A2, b2, ans2, 2, "Matrix1D GE");

    constexpr size_t N3 = 4;
    Matrix1D<N3> A3;
    A3.init_zero();

    // -10   0   1   0 | 1
    //   0   1   0  -1 | 2
    //   3   0  -4  15 | 3     // cols swap
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
    constexpr std::array<double, N3> ans3 = {
        -29. / 224., 481. / 224., -66. / 224., 33. / 224.};
    run_test_with_ans<double, N3, Matrix1D<N3>>(A3, b3, ans3, 3, "Matrix1D GE");

    constexpr size_t N4 = 6;
    Matrix1D<N4> A4;
    A4.init_zero();

    //   1   2   3   0   0   0 |   1
    // -11 -21  22  31   0   0 |   2
    //  66   4   3   5   2   0 |   3
    //   0   3   0   2   4   5 |   4
    //   0   0   5   3   2   0 |   5
    //   0   0   0   1   1  12 |   6
    A4.set(0, 0, 1);
    A4.set(0, 1, 2);
    A4.set(0, 2, 3);
    A4.set(1, 0, -11);
    A4.set(1, 1, -21);
    A4.set(1, 2, 22);
    A4.set(1, 3, 31);
    A4.set(2, 0, 66);
    A4.set(2, 1, 4);
    A4.set(2, 2, 3);
    A4.set(2, 3, 5);
    A4.set(2, 4, 2);
    A4.set(3, 1, 3);
    A4.set(3, 3, 2);
    A4.set(3, 4, 4);
    A4.set(3, 5, 5);
    A4.set(4, 2, 5);
    A4.set(4, 3, 3);
    A4.set(4, 4, 2);
    A4.set(5, 3, 1);
    A4.set(5, 4, 1);
    A4.set(5, 5, 12);
    constexpr std::array<double, N4> b4 = {1, 2, 3, 4, 5, 6};
    constexpr std::array<double, N4> ans4 = {38572. / 342715.,  -77442. / 68543.,
                                             359521. / 342715., -96330. / 68543.,
                                             136092. / 68543.,  30958. / 68543.};
    run_test_with_ans<double, N4, Matrix1D<N4>>(A4, b4, ans4, 4, "Matrix1D GE");

    constexpr size_t N5 = 6;
    MatrixGSL<N5> A5;
    A5.init_zero();

    //   1   2   3   0   0   0 |   1
    // -11 -21  22  31   0   0 |   2
    //  66   4   3   5   2   0 |   3
    //   0   3   0   2   4   5 |   4
    //   0   0   5   3   2   0 |   5
    //   0   0   0   1   1  12 |   6
    A5.set(0, 0, 1);
    A5.set(0, 1, 2);
    A5.set(0, 2, 3);
    A5.set(1, 0, -11);
    A5.set(1, 1, -21);
    A5.set(1, 2, 22);
    A5.set(1, 3, 31);
    A5.set(2, 0, 66);
    A5.set(2, 1, 4);
    A5.set(2, 2, 3);
    A5.set(2, 3, 5);
    A5.set(2, 4, 2);
    A5.set(3, 1, 3);
    A5.set(3, 3, 2);
    A5.set(3, 4, 4);
    A5.set(3, 5, 5);
    A5.set(4, 2, 5);
    A5.set(4, 3, 3);
    A5.set(4, 4, 2);
    A5.set(5, 3, 1);
    A5.set(5, 4, 1);
    A5.set(5, 5, 12);
    constexpr std::array<double, N5> b5 = {1, 2, 3, 4, 5, 6};
    constexpr std::array<double, N5> ans5 = {38572. / 342715.,  -77442. / 68543.,
                                             359521. / 342715., -96330. / 68543.,
                                             136092. / 68543.,  30958. / 68543.};
    run_test_with_ans<double, N5, MatrixGSL<N5>>(A5, b5, ans5, 5, "MatrixGSL GE");

    constexpr size_t N6 = 500;
    Matrix1D<N6> A6;
    std::array<double, N6> b6{};
    A6.init_normal();
    generate_normal<double, N6>(b6, 0.0, 1.0);
    run_test_no_ans(A6, b6, 6, "Matrix1D GE (Norm)", true);

    constexpr size_t N7 = 10;
    MatrixGSL<N7> A7;
    std::array<float, N7> b7{};
    A7.init_normal();
    generate_normal<float, N7>(b7, 0.0, 1.0);

//    run_test_no_ans(A7, b7, 7, "MatrixGSL GE (Norm)", true);


    constexpr size_t n_tests = 12;
    constexpr std::array<size_t, n_tests> Ns {
        50, 100, 150, 200, 250, 300, 350, 400, 450, 500,
        550, 600,
//        650,700, 750, 800, 850, 900, 950, 1000,
//        1050, 1100, 1150, 1200, 1250, 1300, 1350, 1400, 1450, 1500,
//        2200, 2400, 2600, 2800, 3000
    };

    std::array<double, n_tests> norms_ge{};
    norms_ge[ 0] = make_normal_test<double, Ns[ 0], Matrix1D<Ns[ 0]>>(1);
    norms_ge[ 1] = make_normal_test<double, Ns[ 1], Matrix1D<Ns[ 1]>>(2);
    norms_ge[ 2] = make_normal_test<double, Ns[ 2], Matrix1D<Ns[ 2]>>(3);
    norms_ge[ 3] = make_normal_test<double, Ns[ 3], Matrix1D<Ns[ 3]>>(4);
    norms_ge[ 4] = make_normal_test<double, Ns[ 4], Matrix1D<Ns[ 4]>>(5);
    norms_ge[ 5] = make_normal_test<double, Ns[ 5], Matrix1D<Ns[ 5]>>(6);
    norms_ge[ 6] = make_normal_test<double, Ns[ 6], Matrix1D<Ns[ 6]>>(7);
    norms_ge[ 7] = make_normal_test<double, Ns[ 7], Matrix1D<Ns[ 7]>>(8);
    norms_ge[ 8] = make_normal_test<double, Ns[ 8], Matrix1D<Ns[ 8]>>(9);
    norms_ge[ 9] = make_normal_test<double, Ns[ 9], Matrix1D<Ns[ 9]>>(10);
    norms_ge[10] = make_normal_test<double, Ns[10], Matrix1D<Ns[10]>>(11);
    norms_ge[11] = make_normal_test<double, Ns[11], Matrix1D<Ns[11]>>(12);

    std::array<float, n_tests> norms_gsl{};
    norms_gsl[ 0] = make_normal_test<double, Ns[ 0], MatrixGSL<Ns[ 0]>>(1);
    norms_gsl[ 1] = make_normal_test<double, Ns[ 1], MatrixGSL<Ns[ 1]>>(2);
    norms_gsl[ 2] = make_normal_test<double, Ns[ 2], MatrixGSL<Ns[ 2]>>(3);
    norms_gsl[ 3] = make_normal_test<double, Ns[ 3], MatrixGSL<Ns[ 3]>>(4);
    norms_gsl[ 4] = make_normal_test<double, Ns[ 4], MatrixGSL<Ns[ 4]>>(5);
    norms_gsl[ 5] = make_normal_test<double, Ns[ 5], MatrixGSL<Ns[ 5]>>(6);
    norms_gsl[ 6] = make_normal_test<double, Ns[ 6], MatrixGSL<Ns[ 6]>>(7);
    norms_gsl[ 7] = make_normal_test<double, Ns[ 7], MatrixGSL<Ns[ 7]>>(8);
    norms_gsl[ 8] = make_normal_test<double, Ns[ 8], MatrixGSL<Ns[ 8]>>(9);
    norms_gsl[ 9] = make_normal_test<double, Ns[ 9], MatrixGSL<Ns[ 9]>>(10);
    norms_gsl[10] = make_normal_test<double, Ns[10], MatrixGSL<Ns[10]>>(11);
    norms_gsl[11] = make_normal_test<double, Ns[11], MatrixGSL<Ns[11]>>(12);
//    norms[12] = make_normal_test<double, Ns[12], Matrix1D<Ns[12]>>(13);
//    norms[13] = make_normal_test<double, Ns[13], Matrix1D<Ns[13]>>(14);
//    norms[14] = make_normal_test<double, Ns[14], Matrix1D<Ns[14]>>(15);
//    norms[15] = make_normal_test<double, Ns[15], Matrix1D<Ns[15]>>(16);
//    norms[16] = make_normal_test<double, Ns[16], Matrix1D<Ns[16]>>(17);
//    norms[17] = make_normal_test<double, Ns[17], Matrix1D<Ns[17]>>(18);
//    norms[18] = make_normal_test<double, Ns[18], Matrix1D<Ns[18]>>(19);
//    norms[19] = make_normal_test<double, Ns[19], Matrix1D<Ns[19]>>(20);
//    norms[20] = make_normal_test<double, Ns[20], Matrix1D<Ns[20]>>(21);
//    norms[21] = make_normal_test<double, Ns[21], Matrix1D<Ns[21]>>(22);
//    norms[22] = make_normal_test<double, Ns[22], Matrix1D<Ns[22]>>(23);
//    norms[23] = make_normal_test<double, Ns[23], Matrix1D<Ns[23]>>(24);
//    norms[24] = make_normal_test<double, Ns[24], Matrix1D<Ns[24]>>(25);
//    norms[25] = make_normal_test<double, Ns[25], Matrix1D<Ns[25]>>(26);
//    norms[26] = make_normal_test<double, Ns[26], Matrix1D<Ns[26]>>(27);
//    norms[27] = make_normal_test<double, Ns[27], Matrix1D<Ns[27]>>(28);
//    norms[28] = make_normal_test<double, Ns[28], Matrix1D<Ns[28]>>(29);
//    norms[29] = make_normal_test<double, Ns[29], Matrix1D<Ns[29]>>(30);
//    norms[30] = make_normal_test<double, Ns[30], Matrix1D<Ns[30]>>(31);
//    norms[31] = make_normal_test<double, Ns[31], Matrix1D<Ns[31]>>(32);
//    norms[32] = make_normal_test<double, Ns[32], Matrix1D<Ns[32]>>(33);
//    norms[33] = make_normal_test<double, Ns[33], Matrix1D<Ns[33]>>(34);
//    norms[34] = make_normal_test<double, Ns[34], Matrix1D<Ns[34]>>(35);

    save(Ns, norms_ge, norms_gsl, "../test_results/n_vs_norm.csv");
    std::cout << "==================== END TEST SLAE ====================\n\n";
}

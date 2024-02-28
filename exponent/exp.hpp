#pragma once

#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_sf_trig.h>
#include <cstdint>
#include <cstring>
#include <vector>
#include "constants.hpp"
#include "poly.hpp"

namespace ADAAI {

enum class Method { Taylor, Pade, Chebyshev };

template <typename F>
void solveChebyshev(const int N, std::vector<F> &res) {
    double coef[(N + 1) * (N + 1)];
    memset(coef, 0., sizeof(coef));
    for (int k = 0; k < N; k++) {
        coef[k * (N + 1) + k] = -1;
        for (int n = k + 1; n < N + 1; n++) {
            if (n % 2 == 0) {
                if (k % 2 == 1) {
                    coef[k * (N + 1) + n] = 2 * n;
                }
            } else {
                if (k == 0) {
                    coef[k * (N + 1) + n] = n;
                } else if (k % 2 == 0) {
                    coef[k * (N + 1) + n] = 2 * n;
                }
            }
        }
    }
    for (int n = 0; n < N + 1; n++) {
        if (n % 4 == 0) {
            coef[N * (N + 1) + n] = 1;
        } else if (n % 2 == 0) {
            coef[N * (N + 1) + n] = -1;
        }
    }

    double b[N + 1];
    memset(b, 0., sizeof(b));
    b[N] = 1;

    gsl_matrix_view lhs;
    gsl_vector_view rhs;

    lhs = gsl_matrix_view_array(coef, N + 1, N + 1);
    rhs = gsl_vector_view_array(b, N + 1);
    gsl_vector *result = gsl_vector_alloc(N + 1);

    int signum;

    gsl_permutation *lhsPermutation = gsl_permutation_alloc(N + 1);
    gsl_linalg_LU_decomp(&lhs.matrix, lhsPermutation, &signum);

    gsl_linalg_LU_solve(&lhs.matrix, lhsPermutation, &rhs.vector, result);

    for (int i = 0; i <= N; ++i) {
        res[i] = static_cast<F>(gsl_vector_get(result, i));
    }

    gsl_vector_free(result);
    gsl_permutation_free(lhsPermutation);
}

template <typename F, size_t Capacity>
constexpr std::pair<Poly<F, Capacity>, Poly<F, Capacity>>
solvePade(Poly<F, Capacity> T, Poly<F, Capacity> E, size_t n) {
    Poly<F, Capacity> A(std::array<F, Capacity>({1}));
    Poly<F, Capacity> B(std::array<F, Capacity>({0}));
    Poly<F, Capacity> C(std::array<F, Capacity>({0}));
    Poly<F, Capacity> D(std::array<F, Capacity>({1}));

    while (D.deg < n) {
        Poly<F, Capacity> div = ((A * E) + (B * T)) / ((C * E) + (D * T));
        Poly<F, Capacity> C_next = A - div * C;
        Poly<F, Capacity> D_next = B - div * D;
        A = C;
        B = D;
        C = C_next;
        D = D_next;
    }

    return std::make_pair(C * E + D * T, D);
}

// constexpr -- compile-time evaluation
template <Method M = Method::Pade, typename F, size_t Capacity>
constexpr F Exp(F a_x) {
    // F must be floating-point number
    static_assert(std::is_floating_point_v<F>);

    // checking extreme values
    if (a_x > static_cast<F>((INT32_MAX)))
        return std::numeric_limits<F>::infinity();  // almost infinity
    else if (a_x < static_cast<F>((INT32_MIN)))
        return 0.0;  // almost zero

    // simplify calculations by partitioning x into 2 parts
    // idea: represent exp(y) as 2^{n + y0}, where n is an integer
    // and y0 <= 1/2. This form is much more convenient for computation.

    F y = a_x / Ln2<F>();  // y := x / ln2
    F n = NAN;
    F y0 = std::modf(y, &n);  // y := n + y_0

    // make |y0| <= 1/2 preserving n + y0 = y
    if (y0 > 0.5)
        n += 1, y0 -= 1;
    else if (y0 < -0.5)
        n -= 1, y0 += 1;

    // exp(x) = 2^n * 2^y0 = 2^n * exp(y0 * ln2)
    // now, we need to compute y1 := 2^y0 = exp(y0 * ln2) and further use
    // std::ldexp method to get 2^n * y1 = std::ldexp(y1, n).

    // |y0| < 1/2 => |y0 * ln2| < ln2 / 2 => y1 < sqrt(2)
    // exp(x) = \sum_{k=0}^{N} x^k / k! + R_N(x), where holds
    // R_N(x) <= e^ksi * x^{N+1} / (N+1)! and ksi in (0, x)

    // Thus,   R_N(y0 * ln2) <
    //         < e^ksi * (y0 * ln2)^{N+1} / (N+1)! <
    //         < sqrt(2) * (y0 * ln2)^{N+1} / (N+1)!
    // we can choose N so that R_N(y0 * ln2) < atol (absolute tolerance)

    F arg = y0 * Ln2<F>();      // avoid calculating argument several times
    std::vector<F> st = {1.0};  // summing terms for Taylor formula
    F y1 = 0;
    constexpr int N =
        MKExpTaylorOrder<F>();  // compile-time evaluation of order

    if constexpr (M == Method::Taylor) {
        for (int k = 1; k <= N; k++) {
            st.push_back(st.back() * arg / k);
        }
        // compute y1 = 2^y0 = exp(y0 * ln2) via Taylor formula
        // we summarize the terms of the Taylor formula in ascending order of
        // the modulus of the terms to minimize the error of calculations in
        // floating point numbers arithmetic

        for (int i = 0; i < st.size(); ++i)
            y1 += st[st.size() - i - 1];
    } else if constexpr (M == Method::Pade) {
        std::array<F, Capacity> TaylorCoef = {1.0};
        for (int k = 1; k <= N; k++) {
            TaylorCoef[k] = TaylorCoef[k - 1] / k;
        }

        std::array<F, Capacity> monomial = {};
        monomial[N + 1] = 1;

        auto [P, Q] = solvePade(
            Poly<F, Capacity>(TaylorCoef), Poly<F, Capacity>(monomial), N / 2
        );

        y1 = P.eval(arg) / Q.eval(arg);
    } else if constexpr (M == Method::Chebyshev) {
        std::vector<F> c(N + 2);
        solveChebyshev(N + 1, c);
        for (int i = 0; i <= N + 1; ++i) {
            F acos = gsl_complex_arccos_real(arg).dat[0];
            y1 += c[i] * gsl_sf_cos(i * acos);
        }
    }

    return std::ldexp(y1, n);  // return 2^n * y1
}

}  // namespace ADAAI

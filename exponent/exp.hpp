#pragma once

#include <cstdint>
#include <vector>
#include "constants.hpp"
#include "poly.hpp"

namespace ADAAI {

enum class MethodE { Taylor, Pade };

template <typename F, size_t Capacity>
constexpr Poly<F, Capacity>
solvePade(Poly<F, Capacity> T, Poly<F, Capacity> E, size_t n) {
    Poly<F, Capacity> A(std::array<F, Capacity>(1));
    Poly<F, Capacity> B(std::array<F, Capacity>(0));
    Poly<F, Capacity> C(std::array<F, Capacity>(0));
    Poly<F, Capacity> D(std::array<F, Capacity>(1));

    while (D.deg < n) {
        Poly<F, Capacity> div = (A * E + B * T) / (C * E + D * T);
        Poly<F, Capacity> C_next = A + div * C;
        Poly<F, Capacity> D_next = B + div * D;

        A = C;
        B = D;
        C = C_next;
        D = D_next;
    }

    return std::make_pair(C * E + D * T, D);
}

template <typename F>
constexpr F PadeExp(F a_x) {
    Poly<F, PadeNum<F>.size()> P(PadeNum<F>, PadeNum<F>.size());
    Poly<F, PadeDen<F>.size()> Q(PadeDen<F>, PadeDen<F>.size());

    return P.eval(a_x) / Q.eval(a_x);
}

// constexpr -- compile-time evaluation
template <MethodE M = MethodE::Pade, typename F, size_t Capacity>
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

    if constexpr (M == MethodE::Taylor) {
        for (int k = 1; k <= N; k++) {
            st.push_back(st.back() * arg / k);
        }
        // compute y1 = 2^y0 = exp(y0 * ln2) via Taylor formula
        // we summarize the terms of the Taylor formula in ascending order of
        // the modulus of the terms to minimize the error of calculations in
        // floating point numbers arithmetic

        for (int i = 0; i < st.size(); ++i)
            y1 += st[st.size() - i - 1];
    }

    if constexpr (M == MethodE::Pade) {
        std::array<F, Capacity> TaylorCoef = {1.0};
        for (int k = 1; k <= N; k++) {
            TaylorCoef[k] = TaylorCoef[k - 1] * arg / k;
        }
        Poly<F, Capacity> TaylorS(TaylorCoef, N);

        std::array<F, Capacity> monomial = {};
        monomial[N + 1] = 1;

        y1 = PadeExp<F>(arg);
    }

    return std::ldexp(y1, n);  // return 2^n * y1
}

}  // namespace ADAAI

#pragma once

#include <array>
#include <cmath>

namespace ADAAI {

#define M_SQRT2f 1.41421356237309504880f                 /* sqrt(2) */
#define M_SQRT2l 1.414213562373095048801688724209698079L /* sqrt(2) */
#define M_LOG2Ef 1.4426950408889634074f                  /* log_2 e */
#define M_LOG2El 1.442695040888963407359924681001892137L /* log_2 e */

template <typename F>
constexpr inline F Ln2() {
    static_assert(std::is_floating_point_v<F>);
    if (std::is_same<F, float>::value) {
        return 1 / M_LOG2Ef;
    }
    if (std::is_same<F, double>::value) {
        return 1 / M_LOG2E;
    }
    if (std::is_same<F, long double>::value) {
        return 1 / M_LOG2El;
    }
}

template <typename F>
constexpr inline F Sqrt2() {
    static_assert(std::is_floating_point_v<F>);
    if (std::is_same<F, float>::value) {
        return M_SQRT2f;
    }
    if (std::is_same<F, double>::value) {
        return M_SQRT2;
    }
    if (std::is_same<F, long double>::value) {
        return M_SQRT2l;
    }
}

template <typename F>
constexpr inline F Eps = std::numeric_limits<F>::epsilon();

template <typename F>
constexpr inline int MKExpTaylorOrder() {
    F a_atol = 10.0 * Eps<F>;
    F arg = Ln2<F>() / 2;
    F rem = Sqrt2<F>() * arg;
    int k = 1;
    for (; std::abs(rem) > a_atol; k++) {
        rem *= arg / (k + 1);
    }
    return k - 1;
}

}  // namespace ADAAI

#pragma once

#include <array>
#include <cmath>

namespace ADAAI {

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
constexpr inline std::array<F, 2> PadeNum = {1.0, 0.07692307692187728};

template <typename F>
constexpr inline std::array<F, 13> PadeDen = {
    1.0,
    -0.9230769230781227,
    0.4230769230781227,
    -0.12820512820572802,
    0.028846153846353764,
    -0.005128205128255092,
    0.0007478632478732329,
    -9.15750915767523e-05,
    9.53907203930801e-06,
    -8.479175146132913e-07,
    6.359381359694395e-08,
    -3.854170521123118e-09,
    1.6059043838734714e-10};

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

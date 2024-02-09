#pragma once

#include <array>
#include <cmath>

namespace ADAAI {

template <typename F>
constexpr inline F Ln2() {
    return std::log((F)2);
}

template <typename F>
constexpr inline F Sqrt2() {
    return std::sqrt((F)2);
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

}  // namespace ADAAI

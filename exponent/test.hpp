#pragma once

#include <cmath>
#include <functional>

namespace ADAAI {
template <typename T>
// calculates max errors over given interval for a_exponent
std::pair<T, T>
makeTests(T a_l, T a_r, T a_step, std::function<T(T)> a_exponent) {
    T currentX = a_l;
    T absError = 0.0;
    T relError = 0.0;
    while (currentX <= a_r) {
        T stdExp = std::exp(currentX);
        T impExp = a_exponent(currentX);
        T diff = std::abs(stdExp - impExp);
        absError = std::max(diff, absError);
        relError = std::max(diff / stdExp, relError);
        currentX += a_step;
    }
    return std::make_pair(absError, relError);
}
}  // namespace ADAAI

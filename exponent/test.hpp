#pragma once

#include <cmath>
#include <functional>
#include "exp.hpp"

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
        if (currentX < 0) {
            absError = std::max(diff, absError);
        } else {
            relError = std::max(diff / stdExp, relError);
        }
        currentX += a_step;
    }
    return std::make_pair(absError, relError);
}

template <typename F>
// wrapper for tests verbose
void runTests(F a_l, F a_r, F a_step, const std::string &a_typename) {
    const size_t Capacity = 50;

    std::cout << "=== RESULTS FOR " << a_typename << " ===" << std::endl;
    std::cout << "=> EPS is set to be "
              << 10.0 * ADAAI::Eps<float> << std::endl;
    {
        auto [absError, relError] =
            makeTests<F>(a_l, a_r, a_step, Exp<Method::Pade, F, Capacity>);
        std::cout << "=> TAYLOR     | Max absolute error: " << absError
                  << std::endl;
        std::cout << "=> TAYLOR     | Max relative error: " << relError
                  << std::endl;
    }
    {
        auto [absError, relError] =
            makeTests<F>(a_l, a_r, a_step, Exp<Method::Pade, F, Capacity>);
        std::cout << "=> PADE       | Max absolute error: " << absError
                  << std::endl;
        std::cout << "=> PADE       | Max relative error: " << relError
                  << std::endl;
    }
    {
        auto [absError, relError] =
            makeTests<F>(a_l, a_r, a_step, Exp<Method::Chebyshev, F, Capacity>);
        std::cout << "=> CHEBYSHEV  | Max absolute error: " << absError
                  << std::endl;
        std::cout << "=> CHEBYSHEV  | Max relative error: " << relError
                  << std::endl;
    }
    std::cout << std::endl;
}
}  // namespace ADAAI

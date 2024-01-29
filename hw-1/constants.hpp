#ifndef MAIN_CPP_CONSTANTS_HPP
#define MAIN_CPP_CONSTANTS_HPP

#include <cmath>
#include <cfloat>

template<typename F>
constexpr inline F Ln2() {
    return std::log(2);
}

template<typename F>
constexpr inline F Sqrt2() {
    return std::sqrt(2);
}

template<typename F>
constexpr inline F Eps = std::numeric_limits<F>::min();

#endif

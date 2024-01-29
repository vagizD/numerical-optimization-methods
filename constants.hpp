
#ifndef MAIN_CPP_CONSTANTS_HPP
#define MAIN_CPP_CONSTANTS_HPP

#include <cmath>
#include <cfloat>

template<typename F>
constexpr inline F ln2() {                 // natural log of 2
    return std::log(2);
}

template<typename F>
constexpr inline F sqrt2() {
    return std::sqrt(2);            // root of 2
}

template<typename F>
constexpr inline F Eps = std::numeric_limits<F>::epsilon();   // closest value to zero

#endif //MAIN_CPP_CONSTANTS_HPP

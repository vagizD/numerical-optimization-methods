#pragma once
#include <cmath>
#include <vector>


namespace ADAAI {

constexpr inline double INITIAL_STEP  = 1e5;
constexpr inline double INTERVAL_SIZE = 1e-20;

template <typename F>
constexpr inline F Eps = std::numeric_limits<F>::epsilon();

}  // namespace ADAAI

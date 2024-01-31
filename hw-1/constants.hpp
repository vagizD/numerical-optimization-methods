#pragma once

#include <cmath>

namespace ADAAI
{

template <typename F>
constexpr inline F Ln2()
{
  return std::log((F)2);
}

template <typename F>
constexpr inline F Sqrt2()
{
  return std::sqrt((F)2);
}

template <typename F>
constexpr inline F Eps = std::numeric_limits<F>::epsilon();

} // namespace ADAAI

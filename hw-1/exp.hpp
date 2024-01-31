#pragma once

#include "constants.hpp"
#include <cstdint>
#include <vector>

namespace ADAAI
{

// constexpr -- compile-time evaluation
template <typename F>
constexpr F Exp(F a_x, F a_atol = 10.0 * Eps<F>)
{
  // F must be floating-point number
  static_assert(std::is_floating_point_v<F>);

  // checking extreme values
  if (a_x > (F)(INT32_MAX))
    return std::numeric_limits<F>::infinity(); // almost infinity
  else if (a_x < (F)(INT32_MIN))
    return 0.0; // almost zero

  // simplify calculations by partitioning x into 2 parts
  // idea: represent exp(y) as 2^{n + y0}, where n is an integer
  // and y0 <= 1/2. This form is much more convenient for computation.

  F y = a_x / Ln2<F>(); // y := x / ln2
  F n = NAN;
  F y0 = std::modf(y, &n); // y := n + y_0

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

  F arg = y0 * Ln2<F>();     // avoid calculating argument several times
  F rem = Sqrt2<F>() * arg;  // remainder for Taylor formula
  std::vector<F> st = {1.0}; // summing terms for Taylor formula
  for (int k = 1; std::abs(rem) > a_atol; k++)
    {
      rem *= arg / (k + 1);
      st.push_back(st.back() * arg / k);
    }

  // compute y1 = 2^y0 = exp(y0 * ln2) via Taylor formula
  // we summarize the terms of the Taylor formula in ascending order of the
  // modulus of the terms to minimize the error of calculations in floating
  // point numbers arithmetic

  F y1 = 0;
  for (int i = 0; i < st.size(); ++i)
    y1 += st[st.size() - i - 1];
  return std::ldexp(y1, n); // return 2^n * y1
}
} // namespace ADAAI

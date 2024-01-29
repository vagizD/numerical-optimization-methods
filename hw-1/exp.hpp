#include <type_traits>
#include <cmath>
#include <cstdint>

#include "constants.hpp"


template<typename F>                             // F -- should be floating-point number
constexpr F Exp(F x, F atol = 2 * Eps<F>) {      // constexpr -- compile-time evaluation
    static_assert(std::is_floating_point_v<F>);

    // checking extreme values
    if (x > (F)(INT32_MAX))
    {
        return std::numeric_limits<F>::infinity();  // almost infinity
    }
    else if (x < (F)(INT32_MIN))
    {
        return 0.0;  // almost zero
    }

    // simplify calculations by partitioning
    // x into 2 parts
    // Idea: represent exp(y) as 2^{n + y0}, where n is integer
    // and y0 is close to zero. This form is much more
    // convenient for computation.

    //        y := x / ln2 = n + y0
    F y = x / ln2<F>();  // exp(x) = 2^{x / ln2} = 2^n * 2^y0

    F n;                       // fractional part -- |y0| < 1
    F y0 = std::modf(y, &n);  // integral part -- integer

    // make |y0| <= 0.5 preserving n + y0 = y
    if (y0 > 0.5)
    {
        n += 1;
        y0 -= 1;
    }
    else if (y0 < -0.5)
    {
        n -= 1;
        y0 += 1;
    }

    // exp(x) = 2^n * 2^y0 = 2^n * exp(y0 * ln2)

    // now, we need to compute y1 := exp(y0 * ln2) and further use
    // std::ldexp method to get 2^n * y1 = std::ldexp(y1, n).


    // |y0| < 1/2 => |y0 * ln2| < ln2 / 2 => y1 < sqrt(2)
    // exp(X) = \sum_{k=0}^{N} X^k / k! + R_N(X), where holds               (ksi \in (0, X))
    // R_N(X) < e^ksi * X^{N+1} / (N+1)!.

    // Thus,     R_N(y0 * ln2) <
    //         < e^ksi * (y0 * ln2)^{N+1} / (N+1)! <
    //         < sqrt(2) * (y0 * ln2)^{N+1} / (N+1)!
    // we can choose N so that R_N(y0 * ln2) < atol                          (absolute tolerance)

    F arg = y0 * ln2<F>();  // avoid calculating argument several times

    // values for k = 0
    F R = sqrt2<F>() * arg; // remainder
    F st = 1;               // summing term
    F y1 = 1;               // we will compute y1


    for (int k = 1; std::abs(R) > atol; k++) {
        R *= arg / (k+1);
        st *= arg / k;
        y1 += st;
    }

    return std::ldexp(y1, n);
}

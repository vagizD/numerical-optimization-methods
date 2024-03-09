#pragma once

#include <cassert>
#include <cmath>
#include <stdexcept>
#include "aad.h"
#include "enum.h"

template <typename Callable, Derivative D>
double approxStencil3(Callable F, double x, double y, double step = 1e-4) {
    double hx = step, hy = step;
    if (std::abs(x) > 1) {
        hx *= std::abs(x);
    }
    if (std::abs(y) > 1) {
        hy *= std::abs(y);
    }

    switch (D) {
        case Derivative::X:  // Central Difference
            return (F(x + hx, y) - F(x - hx, y)) / (2 * hx);
        case Derivative::Y:  // Central Difference
            return (F(x, y + hy) - F(x, y - hy)) / (2 * hy);
        case Derivative::XX:
            return (F(x + hx, y) - 2 * F(x, y) + F(x - hx, y)) / (hx * hx);
        case Derivative::YY:
            return (F(x, y + hy) - 2 * F(x, y) + F(x, y - hy)) / (hy * hy);
        case Derivative::XY:
            double d_plus = approxStencil3<Callable, Derivative::X>(F, x, y + hy, step);
            double d_minus = approxStencil3<Callable, Derivative::X>(F, x, y - hy, step);
            return (d_plus - d_minus) / (2 * hy);
    }
}

template <typename Callable, Derivative D>
double approxStencil5(Callable F, double x, double y, double step = 1e-4) {
    double hx = step, hy = step;
    if (std::abs(x) > 1) {
        hx *= std::abs(x);
    }
    if (std::abs(y) > 1) {
        hy *= std::abs(y);
    }

    switch (D) {
        case Derivative::X:
            return (-F(x + 2 * hx, y) + 8 * F(x + hx, y) - 8 * F(x - hx, y) +
                    F(x - 2 * hx, y)) /
                   (12 * hx);
        case Derivative::Y:
            return (-F(x, y + 2 * hy) + 8 * F(x, y + hy) - 8 * F(x, y - hy) +
                    F(x, y - 2 * hy)) /
                   (12 * hy);
        case Derivative::XX:
            return (-F(x + 2 * hx, y) + 16 * F(x + hx, y) - 30 * F(x, y) +
                    16 * F(x - hx, y) - F(x - 2 * hx, y)) /
                   (12 * hx * hx);
        case Derivative::YY:
            return (-F(x, y + 2 * hy) + 16 * F(x, y + hy) - 30 * F(x, y) +
                    16 * F(x, y - hy) - F(x, y - 2 * hy)) /
                   (12 * hy * hy);
        case Derivative::XY:
            return (F(x + hx, y + hy) - F(x - hx, y + hy) - F(x + hx, y - hy) +
                    F(x - hx, y - hy)) /
                   (4 * hx * hy);
    }
}

// ======= STENCIL APPROXIMATION METHODS WITH RICHARDSON'S EXTRAPOLATION =======

template <typename Callable, Derivative D, DiffMethod M>
double
approxStencilExtra(Callable F, double x, double y, double step = 1e-4, int n = 10) {
    assert(n % 2 == 0);
    double der_approx, der_approx_grid;
    switch (M) {
        case DiffMethod::Stencil3:
            der_approx = approxStencil3<Callable, D>(F, x, y, step);
            der_approx_grid = approxStencil3<Callable, D>(F, x, y, step / n);
            break;
        case DiffMethod::Stencil5:
            der_approx = approxStencil5<Callable, D>(F, x, y, step);
            der_approx_grid = approxStencil5<Callable, D>(F, x, y, step / n);
            break;
        default:
            throw std::invalid_argument("Wrong DiffMethod provided.");
    }
    return (n * n * der_approx_grid - der_approx) / (n * n - 1);
}

template <Derivative D, DiffMethod M, typename Callable>
double Differentiator(Callable F, double x, double y) {
    if constexpr (M == DiffMethod::Stencil3) {
        return approxStencil3<Callable, D>(F, x, y);
    } else if constexpr (M == DiffMethod::Stencil3Extra) {
        return approxStencilExtra<Callable, D, DiffMethod::Stencil3>(F, x, y);
    } else if constexpr (M == DiffMethod::Stencil5) {
        return approxStencil5<Callable, D>(F, x, y);
    } else if constexpr (M == DiffMethod::Stencil5Extra) {
        return approxStencilExtra<Callable, D, DiffMethod::Stencil5>(F, x, y);
    } else if constexpr (M == DiffMethod::FwdADD) {
        return F(AAD22(Variable::X, x), AAD22(Variable::Y, y)).get_derivative(D);
    }
}

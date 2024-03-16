#pragma once

#include <functional>
#include <iostream>
#include "aad.h"
#include "differentiator.h"
#include "enum.h"

template <Derivative D, DiffMethod M, typename T>
double makeTests(
    std::function<T(T, T)> f,
    std::function<double(double, double)> df,
    double l_x,
    double r_x,
    double step_x,
    double l_y,
    double r_y,
    double step_y
) {
    double max_err = 0;
    for (double x = l_x; x <= r_x; x += step_x) {
        for (double y = l_y; y <= r_y; y += step_y) {
            max_err = std::max(
                max_err,
                std::abs(df(x, y) - Differentiator<D, M, std::function<T(T, T)>>(f, x, y))
            );
        }
    }
    return max_err;
}

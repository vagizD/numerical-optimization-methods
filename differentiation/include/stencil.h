#pragma once

#include "diff.h"

template <typename Callable, Derivative D>
double approxStencil3(Callable F, double x, double y, double step = 1e-4);

template <typename Callable, Derivative D>
double approxStencil5(Callable F, double x, double y, double step = 1e-4);

template <typename Callable, Derivative D, DiffMethod M>
double approxStencilExtra(Callable F,
                          double x,
                          double y,
                          double step = 1e-4,
                          int n = 10);

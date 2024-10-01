#pragma once
#include "slae.h"
#include "system_params.h"
#include "utils.h"

template <typename С, const size_t m, const size_t n>
void option_price_stepper(
    const SystemParams &a_sp,
    С &A,
    std::array<double, (m + 1) * (n + 1)> &x
) {
    size_t k = a_sp.m_p;
    std::array<double, (m + 1) * (n + 1)> b = {};

    init_b<m, n>(x, a_sp);
    init_rlhs<m, n>(A, b, x, a_sp, k);
    while (k >= 1) {
        solve_slae(A, b, x);
        update_b<m, n>(b, x, a_sp, k);
        k--;
    }
}

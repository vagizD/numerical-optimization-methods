#pragma once
#include "container.h"
#include "slae.h"
#include "system_params.h"
#include "utils.h"

template <typename Cont, ContainerType type, const size_t m, const size_t n>
void option_price_stepper(
    const SystemParams &a_sp,
    Cont &A,
    std::array<double, (m + 1) * (n + 1)> &x
) {
    // A.init_zero();
    size_t k = a_sp.m_p;
    std::array<double, (m + 1) * (n + 1)> b = {};

    init_b<m, n>(b, a_sp);
    init_rlhs<m, n>(A, b, x, a_sp, k);
    // update_b<m, n>(b, x, a_sp, k);
    while (k >= 1) {
        if constexpr (type == ContainerType::MatrixSparse) {
            solve_slae_sparse(A, b, x);
        } else {
            solve_slae(A, b, x);
        }
        // k--;
        update_b<m, n>(b, x, a_sp, k);
        k--;
    }
}

#pragma once
#include <array>
#include <algorithm>
#include <cassert>


template <typename F, const size_t N>
F vector_norm(std::array<F, N> &v) {
    assert(N > 0);

    F max_element = std::abs(v[0]);
    for (size_t i = 1; i < N; ++i) {
        max_element = std::max(max_element, std::abs(v[i]));
    }
    return max_element;
}


template <typename F, typename C, const size_t N>
void sqmatrix_on_vector(const C& A, const std::array<F, N> &v, std::array<F, N> &res) {
    for (size_t i = 0; i < N; ++i) {
        res[i] = 0;
        for (size_t j = 0; j < N; ++j) {
            res[i] += A.get(i, j) * v[j];
        }
    }
}


template <typename F, const size_t N>
void vectors_diff(
    const std::array<F, N> &v,
    const std::array<F, N> &u,
    std::array<F, N> &res
) {
    for (size_t i = 0; i < N; ++i) {
        res[i] = v[i] - u[i];
    }
}

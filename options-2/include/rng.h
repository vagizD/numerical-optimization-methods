#pragma once
#include <random>

template <typename F, const size_t N>
void generate_normal(std::array<F, N> &v, F mean, F stddev) {
    std::random_device rd{};
    std::mt19937 gen{rd()};
    std::normal_distribution d(mean, stddev);

    for (size_t i = 0; i < N; ++i) {
        v[i] = d(gen);
    }
}

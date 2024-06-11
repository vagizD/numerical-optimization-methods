#pragma once

#include <array>
#include <cstdlib>
#include "constants.h"

template <size_t N>
void backtrack(
    std::array<double, N> &A,
    std::array<double, N> &B,
    std::array<double, N> &D,
    std::array<double, N> &F,
    std::array<double, N> &c_res
) {
    for (size_t i = 2; i <= N - 1; ++i) {
        B[i] = B[i] - A[i] * D[i - 1] / B[i - 1];
    }

    c_res[N - 2] = F[N - 2] / B[N - 2];
    for (size_t i = N - 3; i > 0; --i) {
        c_res[i] = (F[i] - D[i] * c_res[i + 1]) / B[i];
    }
}

template <size_t N, size_t M>
void make_step(
    double cur_tau,
    std::array<double, N> &c_cur,
    std::array<double, N> &c_res
) {
    std::array<double, N> A;
    std::array<double, N> B;
    std::array<double, N> D;
    std::array<double, N> F;

    double h = 1.0 / N;
    double d_tau = 1.0 / M;
    for (size_t i = 0; i < N; ++i) {
        double si = (double)i * S_max / (double)N;
        double a = risk(cur_tau) * si;
        double b = volatility(cur_tau) * si;
        double g = -risk(cur_tau);
        A[i] = (b * b) / (2.0 * h * h) - a / h;
        B[i] = g - (b * b) / (h * h) - 1.0 / d_tau;
        D[i] = (b * b) / (2.0 * h * h) + a / h;
        F[i] = -c_cur[i] / d_tau;
    }

    // F[1] = F[1] - A[1] * c_cur[0];                  // c_cur[0] == 0
    F[N - 2] =
        F[N - 2] -
        D[N - 2] *
            c_cur[N - 1];  // c_cur[N - 1] == (S_max - (K * std::exp(-integrate_risk(tau))

    backtrack(A, B, D, F, c_res);
}

template <size_t N, size_t M>
void go_implicit(
    const std::array<double, N> &c_init,
    const std::array<double, M> &tau,
    std::array<double, N> &c_res
) {
    std::array<double, N> c_cur = c_init;
    for (size_t j = 0; j < M; ++j) {
        c_cur[0] = 0;
        c_cur[N - 1] = S_max - (double)K * std::exp(-integrate_risk(tau[j]));

        make_step<N, M>(tau[j], c_cur, c_res);
        c_cur = c_res;
    }
    c_res[0] = 0;
    c_res[N - 1] = S_max - (double)K * std::exp(-integrate_risk(tau[M - 1]));
}

#pragma once
#include <array>
#include <cmath>
#include "container.h"
#include "system_params.h"

// ===================== SLAE: Ax = b ============================

// map 2D index to 1D
template <const size_t n>
size_t map2to1(const size_t i, const size_t j) {
    return (n + 1) * i + j;
}

// S is discretized with index j = 0...n
template <const size_t n>
double S(const size_t j, const double S_max) {
    return static_cast<double>(j) / static_cast<double>(n) * S_max;
}

// V is discretized with index i = 0...m
template <const size_t m>
double V(const size_t i, const double V_max) {
    return static_cast<double>(i) / static_cast<double>(m) * V_max;
}

template <const size_t n>
void fill_index_matrix(
    const size_t i,
    const size_t j,
    std::array<std::array<size_t, 3>, 3> &a_im
) {
    for (int dx = -1; dx <= 1; ++dx) {
        for (int dy = -1; dy <= 1; ++dy) {
            a_im[1 + dx][1 + dy] = map2to1<n>(i + dx, j + dy);
        }
    }
}

// Function to:
//   * init A and init b
//   * init A and update b
template <const size_t m, const size_t n>
void init_rlhs(
    Container &a_A,
    std::array<double, (n + 1) * (m + 1)> &a_b,
    std::array<double, (n + 1) * (m + 1)> &a_b_prev,
    const SystemParams &a_sp,
    const size_t k
) {
    a_A.init_zero();
    std::fill(a_b.begin(), a_b.end(), 0);

    // INNER PART EQUATION
    // C(Sj, Vi, t{k-1}) - tau * L * C(Sj, Vi, t{k-1}) = C(Sj, Vi, tk)
    // (1)               - tau * L * (2)               = (3)

    const double tau = a_sp.m_T / static_cast<double>(a_sp.m_p);
    const double mu = a_sp.m_mp.m_r - a_sp.m_mp.m_d;
    const double time = static_cast<double>(a_sp.m_p - k) * tau;

    // LEFT BOUNDARY
    // j = 0 -> S = 0 -> C(0, V, t) = 0
    for (size_t i = 0; i < m + 1; ++i) {
        size_t i0 = map2to1<n>(i, 0);
        a_A.add(i0, i0, 1);
        a_b[i0] = 0;
    }

    // LOW BOUNDARY ((0, 0) in LEFT)
    // i = 0 -> V = 0 -> C(S, 0, t) = 0 if S < K* else ...
    for (size_t j = 1; j < n + 1; ++j) {
        size_t _0j = map2to1<n>(0, j);
        const double Sj = S<n>(j, a_sp.m_S_max);
        const double diff = std::exp(-a_sp.m_mp.m_d * time) * Sj -
                            std::exp(-a_sp.m_mp.m_r * time) * a_sp.m_K;
        const double K_star = std::exp(-mu * time) * a_sp.m_K;

        a_A.add(_0j, _0j, 1);
        if (Sj < K_star) {
            a_b[_0j] = 0;
        } else {
            a_b[_0j] = diff;
        }
    }

    // RIGHT BOUNDARY ((0, n) in LOW)
    // j = n -> S = S_max -> C(S_max, V, t) = S* - K*
    for (size_t i = 1; i < m + 1; ++i) {
        size_t in = map2to1<n>(i, n);
        const double diff = std::exp(-a_sp.m_mp.m_d * time) * a_sp.m_S_max -
                            std::exp(-a_sp.m_mp.m_r * time) * a_sp.m_K;

        a_A.add(in, in, 1);
        a_b[in] = diff;
    }

    // UPPER BOUNDARY ((m, n) in RIGHT, (m, 0) in LEFT)
    // i = m -> V = V_max
    // 0 ... 1 ... 0 0 = 0 -> x_i = 0
    // 0 ... 0 ... 1 1 = 0 -> x_m = x_{m-1}
    for (size_t j = 1; j < n; ++j) {
        const size_t mj = map2to1<n>(m, j);
        const size_t mmj = map2to1<n>(m - 1, j);

        a_A.add(mj, mj, 1);
        a_A.add(mj, mmj, 1);
        a_b[mj] = 0;
    }

    // INNER PART
    for (size_t i = 1; i < m; ++i) {
        for (size_t j = 1; j < n; ++j) {
            std::array<std::array<size_t, 3>, 3> im{};
            fill_index_matrix<n>(i, j, im);

            const size_t ij = map2to1<n>(i, j);
            const double Sj = S<n>(j, a_sp.m_S_max);
            const double Vi = V<m>(i, a_sp.m_V_max);

            a_A.add(ij, im[1][1], 1);  // (1)

            double tmp = mu * Sj / (2 * a_sp.m_hs);
            tmp *= -tau;
            a_A.add(ij, im[1][2], -tmp);  // (2)
            a_A.add(ij, im[1][0], tmp);   // (2)

            tmp = a_sp.m_mp.m_sigma * a_sp.m_mp.m_sigma / 2 *
                  std::pow(Sj, 2 * a_sp.m_mp.m_beta) / (a_sp.m_hs * a_sp.m_hs);
            tmp *= -tau;
            a_A.add(ij, im[1][2], tmp);       // (2)
            a_A.add(ij, im[1][1], -2 * tmp);  // (2)
            a_A.add(ij, im[1][0], tmp);       // (2)

            tmp = a_sp.m_mp.m_kappa * (a_sp.m_mp.m_theta - Vi) / (2 * a_sp.m_hv);
            tmp *= -tau;
            a_A.add(ij, im[2][1], tmp);   // (2)
            a_A.add(ij, im[0][1], -tmp);  // (2)

            tmp = a_sp.m_mp.m_epsilon * a_sp.m_mp.m_epsilon / 2 * Vi /
                  (a_sp.m_hv * a_sp.m_hv);
            tmp *= -tau;
            a_A.add(ij, im[2][1], tmp);       // (2)
            a_A.add(ij, im[1][1], -2 * tmp);  // (2)
            a_A.add(ij, im[0][1], tmp);       // (2)

            tmp = a_sp.m_mp.m_sigma * a_sp.m_mp.m_epsilon * a_sp.m_mp.m_rho *
                  std::pow(Sj, a_sp.m_mp.m_beta) * Vi / (4 * a_sp.m_hs * a_sp.m_hv);
            tmp *= -tau;
            a_A.add(ij, im[2][2], tmp);   // (2)
            a_A.add(ij, im[2][0], -tmp);  // (2)
            a_A.add(ij, im[0][2], -tmp);  // (2)
            a_A.add(ij, im[0][0], tmp);   // (2)

            tmp = a_sp.m_mp.m_r;
            tmp *= -tau;
            a_A.add(ij, im[1][1], -tmp);  // (2)

            // NOT a_A_prev[ij][ij] since it is a coef of Cij,
            // we need previous value of Cij
            a_b[ij] += a_b_prev[ij];  // (3)
        }
    }
}

// Function to init b for t=0
template <const size_t m, const size_t n>
void init_b(std::array<double, (n + 1) * (m + 1)> &a_b, const SystemParams &a_sp) {
    for (size_t j = 0; j < n + 1; ++j) {
        const double Sj = S<n>(j, a_sp.m_S_max);
        for (size_t i = 0; i < m + 1; ++i) {
            size_t ij = map2to1<n>(i, j);
            a_b[ij] = std::max(Sj - a_sp.m_K, 0.0);
        }
    }
}

template <const size_t m, const size_t n>
void update_b(
    std::array<double, (n + 1) * (m + 1)> &a_b,
    std::array<double, (n + 1) * (m + 1)> &a_b_prev,
    const SystemParams &a_sp,
    const size_t k
) {
    std::fill(a_b.begin(), a_b.end(), 0);

    const double tau = a_sp.m_T / static_cast<double>(a_sp.m_p);
    const double mu = a_sp.m_mp.m_r - a_sp.m_mp.m_d;
    const double time = static_cast<double>(a_sp.m_p - k) * tau;

    // LEFT BOUNDARY
    for (size_t i = 0; i < m + 1; ++i) {
        size_t i0 = map2to1<n>(i, 0);
        a_b[i0] = 0;
    }

    // LOW BOUNDARY ((0, 0) in LEFT)
    for (size_t j = 1; j < n + 1; ++j) {
        size_t _0j = map2to1<n>(0, j);
        const double Sj = S<n>(j, a_sp.m_S_max);
        const double diff = std::exp(-a_sp.m_mp.m_d * time) * Sj -
                            std::exp(-a_sp.m_mp.m_r * time) * a_sp.m_K;
        const double K_star = std::exp(-mu * time) * a_sp.m_K;

        if (Sj < K_star) {
            a_b[_0j] = 0;
        } else {
            a_b[_0j] = diff;
        }
    }

    // RIGHT BOUNDARY ((0, n) in LOW)
    for (size_t i = 1; i < m + 1; ++i) {
        size_t in = map2to1<n>(i, n);
        const double diff = std::exp(-a_sp.m_mp.m_d * time) * a_sp.m_S_max -
                            std::exp(-a_sp.m_mp.m_r * time) * a_sp.m_K;

        a_b[in] = diff;
    }

    // UPPER BOUNDARY ((m, n) in RIGHT, (m, 0) in LEFT)
    for (size_t j = 1; j < n; ++j) {
        const size_t mj = map2to1<n>(m, j);
        a_b[mj] = 0;
    }

    // INNER PART
    for (size_t i = 1; i < m; ++i) {
        for (size_t j = 1; j < n; ++j) {
            const size_t ij = map2to1<n>(i, j);
            a_b[ij] += a_b_prev[ij];
        }
    }
}

// PSEUDOCODE

// b_prev = 0;
// init_b(b_prev)
// A, b = 0, 0
// init_rlhs(A, b, b_prev)
// x = 0
// for (k = p, k >= 0):
//     solve_slae(A, b, x) -- solve for x
//     update_b(b, x, k)   -- update b with values of x and boundary conditions
// return x

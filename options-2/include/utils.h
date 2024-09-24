#pragma once
#include <array>
#include <cmath>
#include <tuple>
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
            a_im[i + dx][j + dy] = map2to1<n>(i + dx, j + dy);
        }
    }
}

// // Function to update matrix A and b (A is 2D)
// template<const size_t m, const size_t n>
// void update_rlhs(
//     std::array<std::array<double, (n + 1) * (m + 1)>, (n + 1) * (m + 1)> &a_A,
//     std::array<double, (n + 1) * (m + 1)> &a_b,
//     std::array<double, (n + 1) * (m + 1)> &a_b_prev,
//     const SystemParams &a_sp,
//     const size_t k
//     ) {
//     std::fill(a_A.begin(), a_A.end(), 0);
//     std::fill(a_b.begin(), a_b.end(), 0);
//
//     // INNER PART EQUATION
//     // C(Sj, Vi, t{k-1}) - time * L * C(Sj, Vi, t{k-1}) = C(Sj, Vi, tk)
//     // (1)               - time * L * (2)               = (3)
//
//     const double tau = a_sp.m_T / static_cast<double>(a_sp.m_p);
//     const double mu = a_sp.m_mp.m_r - a_sp.m_mp.m_d;
//     const double time = static_cast<double>(a_sp.m_p - k) * tau;
//
//     // LEFT BOUNDARY
//     // j = 0 -> S = 0 -> C(0, V, t) = 0
//     for (size_t i = 0; i < m + 1; ++i) {
//         size_t i0 = map2to1<n>(i, 0);
//         a_A[i0][i0] += 1;
//         a_b[i0] = 0;
//     }
//
//     // LOW BOUNDARY ((0, 0) in LEFT)
//     // i = 0 -> V = 0 -> C(S, 0, t) = 0 if S < K* else ...
//     for (size_t j = 1; j < n + 1; ++j) {
//         size_t _0j = map2to1<n>(0, j);
//         const double Sj = S<n>(j, a_sp.m_S_max);
//         const double diff = std::exp(- a_sp.m_mp.m_d * time) * Sj - std::exp(-
//         a_sp.m_mp.m_r * time) * a_sp.m_K; const double K_star = std::exp(- mu * time) *
//         a_sp.m_K;
//
//         a_A[_0j][_0j] += 1;
//         if (Sj < K_star) {
//             a_b[_0j] = 0;
//         } else {
//             a_b[_0j] = diff;
//         }
//     }
//
//     // RIGHT BOUNDARY ((0, n) in LOW)
//     // j = n -> S = S_max -> C(S_max, V, t) = S* - K*
//     for (size_t i = 1; i < m + 1; ++i) {
//         size_t in = map2to1<n>(i, n);
//         const double diff = std::exp(- a_sp.m_mp.m_d * time) * a_sp.m_S_max -
//         std::exp(- a_sp.m_mp.m_r * time) * a_sp.m_K;
//
//         a_A[in][in] += 1;
//         a_b[in] = diff;
//     }
//
//     // UPPER BOUNDARY ((m, n) in RIGHT, (m, 0) in LEFT)
//     // i = m -> V = V_max
//     // 0 ... 1 ... 0 0 = 0 -> x_i = 0
//     // 0 ... 0 ... 1 1 = 0 -> x_m = x_{m-1}
//     for (size_t j = 1; j < n; ++j) {
//         size_t mj  = map2to1<n>(m, j);
//         size_t mmj = map2to1<n>(m - 1, j);
//
//         a_A[mj][mj]  += 1;
//         a_A[mj][mmj] += 1;
//         a_b[mj] = 0;
//     }
//
//     // INNER PART
//     for (size_t i = 1; i < m; ++i) {
//         for (size_t j = 1; j < n; ++j) {
//             std::array<std::array<size_t, 3>, 3> im{};
//             fill_index_matrix<n>(i, j, im);
//
//             const size_t ij = map2to1<n>(i, j);
//             const double Sj = S<n>(j, a_sp.m_S_max);
//             const double Vi = V<m>(i, a_sp.m_V_max);
//
//             a_A[ij][ij] += 1;               // (1)
//
//             double tmp = mu * Sj / (2 * a_sp.m_hs);
//             tmp *= -time;
//             a_A[ij][im[1][2]] -= tmp;       // (2)
//             a_A[ij][im[1][0]] += tmp;       // (2)
//
//             tmp = a_sp.m_mp.m_sigma * a_sp.m_mp.m_sigma / 2 * std::pow(Sj, 2 *
//             a_sp.m_mp.m_beta) / (a_sp.m_hs * a_sp.m_hs); tmp *= -time;
//             a_A[ij][im[1][2]] += tmp;       // (2)
//             a_A[ij][im[1][1]] -= 2 * tmp;   // (2)
//             a_A[ij][im[1][0]] += tmp;       // (2)
//
//             tmp = a_sp.m_mp.m_kappa * (a_sp.m_mp.m_theta - Vi) / (2 * a_sp.m_hv);
//             tmp *= -time;
//             a_A[ij][im[2][1]] += tmp;       // (2)
//             a_A[ij][im[0][1]] -= tmp;       // (2)
//
//             tmp = a_sp.m_mp.m_epsilon * a_sp.m_mp.m_epsilon / 2 * Vi / (a_sp.m_hv *
//             a_sp.m_hv); tmp *= -time; a_A[ij][im[2][1]] += tmp;       // (2)
//             a_A[ij][im[1][1]] -= 2 * tmp;   // (2)
//             a_A[ij][im[0][1]] += tmp;       // (2)
//
//             tmp = a_sp.m_mp.m_sigma * a_sp.m_mp.m_epsilon * a_sp.m_mp.m_rho *
//             std::pow(Sj, a_sp.m_mp.m_beta) * Vi / (4 * a_sp.m_hs * a_sp.m_hv); tmp *=
//             -time; a_A[ij][im[2][2]] += tmp;       // (2) a_A[ij][im[2][0]] -= tmp; //
//             (2) a_A[ij][im[0][2]] -= tmp;       // (2) a_A[ij][im[0][0]] += tmp; // (2)
//
//             tmp = a_sp.m_mp.m_r;
//             tmp *= -time;
//             a_A[ij][im[1][1]] -= tmp;       // (2)
//
//             if (k == a_sp.m_p) {
//                 // NOT a_A_prev[ij][ij] since it is a coef of Cij, we need previous
//                 value of Cij a_b[ij] += a_b_prev[ij];    // (3)
//             }
//         }
//     }
// }

// Function to update matrix A and b
template <const size_t m, const size_t n>
void update_rlhs(
    std::array<double, (n + 1) * (m + 1) * (n + 1) * (m + 1)> &a_A,
    std::array<double, (n + 1) * (m + 1)> &a_b,
    std::array<double, (n + 1) * (m + 1)> &a_b_prev,
    const SystemParams &a_sp,
    const size_t k
) {
    std::fill(a_A.begin(), a_A.end(), 0);
    std::fill(a_b.begin(), a_b.end(), 0);

    // INNER PART EQUATION
    // C(Sj, Vi, t{k-1}) - time * L * C(Sj, Vi, t{k-1}) = C(Sj, Vi, tk)
    // (1)               - time * L * (2)               = (3)

    const size_t N = (n + 1) * (m + 1);
    const double tau = a_sp.m_T / static_cast<double>(a_sp.m_p);
    const double mu = a_sp.m_mp.m_r - a_sp.m_mp.m_d;
    const double time = static_cast<double>(a_sp.m_p - k) * tau;

    // LEFT BOUNDARY
    // j = 0 -> S = 0 -> C(0, V, t) = 0
    for (size_t i = 0; i < m + 1; ++i) {
        size_t i0 = map2to1<n>(i, 0);
        a_A[map2to1<N>(i0, i0)] += 1;
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

        a_A[map2to1<N>(_0j, _0j)] += 1;
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

        a_A[map2to1<N>(in, in)] += 1;
        a_b[in] = diff;
    }

    // UPPER BOUNDARY ((m, n) in RIGHT, (m, 0) in LEFT)
    // i = m -> V = V_max
    // 0 ... 1 ... 0 0 = 0 -> x_i = 0
    // 0 ... 0 ... 1 1 = 0 -> x_m = x_{m-1}
    for (size_t j = 1; j < n; ++j) {
        size_t mj = map2to1<n>(m, j);
        size_t mmj = map2to1<n>(m - 1, j);

        a_A[map2to1<N>(mj, mj)] += 1;
        a_A[map2to1<N>(mj, mmj)] += 1;
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

            a_A[map2to1<N>(ij, ij)] += 1;  // (1)

            double tmp = mu * Sj / (2 * a_sp.m_hs);
            tmp *= -time;
            a_A[map2to1<N>(ij, im[1][2])] -= tmp;  // (2)
            a_A[map2to1<N>(ij, im[1][0])] += tmp;  // (2)

            tmp = a_sp.m_mp.m_sigma * a_sp.m_mp.m_sigma / 2 *
                  std::pow(Sj, 2 * a_sp.m_mp.m_beta) / (a_sp.m_hs * a_sp.m_hs);
            tmp *= -time;
            a_A[map2to1<N>(ij, im[1][2])] += tmp;      // (2)
            a_A[map2to1<N>(ij, im[1][1])] -= 2 * tmp;  // (2)
            a_A[map2to1<N>(ij, im[1][0])] += tmp;      // (2)

            tmp = a_sp.m_mp.m_kappa * (a_sp.m_mp.m_theta - Vi) / (2 * a_sp.m_hv);
            tmp *= -time;
            a_A[map2to1<N>(ij, im[2][1])] += tmp;  // (2)
            a_A[map2to1<N>(ij, im[0][1])] -= tmp;  // (2)

            tmp = a_sp.m_mp.m_epsilon * a_sp.m_mp.m_epsilon / 2 * Vi /
                  (a_sp.m_hv * a_sp.m_hv);
            tmp *= -time;
            a_A[map2to1<N>(ij, im[2][1])] += tmp;      // (2)
            a_A[map2to1<N>(ij, im[1][1])] -= 2 * tmp;  // (2)
            a_A[map2to1<N>(ij, im[0][1])] += tmp;      // (2)

            tmp = a_sp.m_mp.m_sigma * a_sp.m_mp.m_epsilon * a_sp.m_mp.m_rho *
                  std::pow(Sj, a_sp.m_mp.m_beta) * Vi / (4 * a_sp.m_hs * a_sp.m_hv);
            tmp *= -time;
            a_A[map2to1<N>(ij, im[2][2])] += tmp;  // (2)
            a_A[map2to1<N>(ij, im[2][0])] -= tmp;  // (2)
            a_A[map2to1<N>(ij, im[0][2])] -= tmp;  // (2)
            a_A[map2to1<N>(ij, im[0][0])] += tmp;  // (2)

            tmp = a_sp.m_mp.m_r;
            tmp *= -time;
            a_A[map2to1<N>(ij, im[1][1])] -= tmp;  // (2)

            if (k == a_sp.m_p) {
                // NOT a_A_prev[ij][ij] since it is a coef of Cij,
                // we need previous value of Cij
                a_b[ij] += a_b_prev[ij];  // (3)
            }
        }
    }
}

// Function to init A and b for k = p
template <const size_t m, const size_t n>
void init_rlhs(
    std::array<std::array<double, (n + 1) * (m + 1)>, (n + 1) * (m + 1)> &a_A,
    std::array<double, (n + 1) * (m + 1)> &a_b,
    const SystemParams &a_sp
) {
    const size_t N = (n + 1) * (m + 1);
    std::array<double, N> b_prev;

    for (size_t j = 0; j < n + 1; ++j) {
        const double Sj = S<n>(j, a_sp.m_S_max);
        for (size_t i = 0; i < m + 1; ++i) {
            size_t ij = map2to1<n>(i, j);
            b_prev[ij] = std::max(Sj - a_sp.m_K, 0.0);
        }
    }

    update_rlhs<m, n>(a_A, a_b, b_prev, a_sp, a_sp.m_p);  // maybe +1 in time = ... ????
}

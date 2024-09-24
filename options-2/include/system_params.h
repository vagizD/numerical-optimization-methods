#pragma once
#include <cassert>

class ModelParams {
public:
    double m_r;        // risk-interest rate
    double m_d;        // dividend rate
    double m_sigma;    // volatility coefficient
    double m_kappa;    // stochastic volatility coefficient
    double m_beta;     // constant elasticity of variance
    double m_epsilon;  // stochastic volatility coefficient
    double m_rho;      // correlation of brownian differentials
    double m_theta;    //

    constexpr ModelParams(
        const double a_r,
        const double a_d,
        const double a_sigma,
        const double a_kappa,
        const double a_beta,
        const double a_epsilon,
        const double a_rho,
        const double a_theta
    )
        : m_r(a_r),
          m_d(a_d),
          m_sigma(a_sigma),
          m_kappa(a_kappa),
          m_beta(a_beta),
          m_epsilon(a_epsilon),
          m_rho(a_rho),
          m_theta(a_theta) {
        assert((0 <= a_rho) && (a_rho < 1));
        assert(a_kappa > 0);
        assert(a_theta > 0);
        assert(a_epsilon > 0);
        assert(a_sigma > 0);
        assert((0 < a_beta) && (a_beta <= 1));
    }
};

class SystemParams {
public:
    // Model Parameters
    const ModelParams m_mp;

    // Task-Specific Parameters
    const double m_T;
    const double m_M;
    const double m_K;
    const double m_S_0;
    const double m_V_0;
    size_t m_S_0_j = 0;
    size_t m_V_0_i = 0;
    double m_S_max = 0;
    double m_V_max = 0;

    // Discretization Parameters
    const size_t m_n;
    const size_t m_m;
    const size_t m_p;
    double m_hs = 0;
    double m_hv = 0;

    SystemParams(
        const double a_r,
        const double a_d,
        const double a_sigma,
        const double a_kappa,
        const double a_beta,
        const double a_epsilon,
        const double a_rho,
        const double a_theta,
        const double a_S_0,
        const double a_V_0,
        const double a_T,
        const double a_M,
        const double a_K,
        const size_t a_n,
        const size_t a_m,
        const size_t a_p
    )
        : m_mp(a_r, a_d, a_sigma, a_kappa, a_beta, a_epsilon, a_rho, a_theta),
          m_T(a_T),
          m_M(a_M),
          m_K(a_K),
          m_S_0(a_S_0),
          m_V_0(a_V_0),
          m_n(a_n),
          m_m(a_m),
          m_p(a_p) {
        auto [j, i, S_max, V_max] =
            adjust_grid_boundary(m_mp, m_n, m_m, m_T, m_M, m_S_0, m_V_0, m_K);
        m_V_0_i = i;
        m_S_0_j = j;
        m_S_max = S_max;
        m_V_max = V_max;
        m_hs = S_max / static_cast<double>(m_n);
        m_hv = V_max / static_cast<double>(m_m);
    }

    // computes S_max, V_max estimates and adjusts them such that S_0 and V_0 lie on a
    // grid S_max and V_max are both increased such that a = S_0 * n / S_max becomes the
    // next closest integer, i.e. i* = [a]
    static std::tuple<std::size_t, std::size_t, double, double> adjust_grid_boundary(
        const ModelParams &a_mp,
        size_t a_n,
        size_t a_m,
        double a_T,
        double a_M,
        double a_S_0,
        double a_V_0,
        double a_K
    );
};

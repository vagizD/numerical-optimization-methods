#pragma once

#include <array>
#include "constants.h"
#include "stepper.h"

double calc(int i, double tau, double c_prev, double c_cur, double c_next) {
    double r = risk(tau);
    double v = volatility(tau);
    double si = (double)i * S_max / (double)N;
    double dc = (c_next - c_prev) / (2.0 * STEP);
    double ddc = (c_next - 2.0 * c_cur + c_prev) / (STEP * STEP);
    return r * si * dc + v * v * si * si * ddc / 2.0 - r * c_cur;
}

template <size_t N_>
class RHS {
public:
    constexpr static int N = N_;
    RHS() = default;

    void differentiate(
        double tau,
        const std::array<double, N> &c_cur,
        std::array<double, N> &dc
    ) const {
        dc = {};
        assert(N >= 2);
        for (int i = 1; i < N - 1; ++i) {
            dc[i] = calc(i, tau, c_cur[i - 1], c_cur[i], c_cur[i + 1]);
        }
    }
};

template <size_t N>
void go(const std::array<double, N> &c_init, std::array<double, N> &c_res) {
    double tau = 0;
    RHS<N> rhs = {};
    TimeStepper_RKF45<RHS<N>> stepper(&rhs);
    std::array<double, N> c_cur = c_init;
    double a_h = STEP;
    while (tau + STEP <= 1.0) {
        double a_t = tau;
        c_cur[0] = 0;
        c_cur[N - 1] = S_max - (double)K * std::exp(-integrate_risk(tau));
        stepper.make_step(a_t, a_h, c_cur, c_res);
        tau += STEP;
        c_cur = c_res;
    }
}

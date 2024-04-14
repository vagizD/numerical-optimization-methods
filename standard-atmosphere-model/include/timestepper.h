#pragma once

#include <array>
#include <utility>
#include <cmath>
#include "constants.h"

template <typename RHS>
class TimeStepper_RKF45 {
private:
    RHS const *const m_rhs;
    constexpr static int N = RHS::N;
    double eps = 1e-6;

public:
    TimeStepper_RKF45(RHS const *a_rhs) : m_rhs(a_rhs){};

    void operator()(
        double& a_t,
        double& a_h,  // suggested step
        std::array<double, RHS::N> &a_y,
        std::array<double, RHS::N> &a_y_next
    ) {
        bool flag = true;
        while (flag) {
            // RKF45 method
            std::array<double, RHS::N> rhs_result;
            // k1
            m_rhs(a_t, a_y, rhs_result);
            std::array<double, RHS::N> k1 = {};
            for (int i = 0; i < RHS::N; i++) {
                k1[i] = a_h * rhs_result[i];
            }
            // k2
            std::array<double, RHS::N> y_for_k2 = {};
            for (int i = 0; i < RHS::N; i++) {
                y_for_k2[i] = a_y[i] + k2_consts[0] * k1[i];
            }
            m_rhs(a_t + 1.0 / 4.0 * a_h, y_for_k2, rhs_result);
            std::array<double, RHS::N> k2 = {};
            for (int i = 0; i < RHS::N; i++) {
                k2[i] = a_h * rhs_result[i];
            }
            // k3
            std::array<double, RHS::N> y_for_k3 = {};
            for (int i = 0; i < RHS::N; i++) {
                y_for_k3[i] = a_y[i] + k3_consts[0] * k1[i] + k3_consts[1] * k2[i];
            }
            m_rhs(a_t + 3.0 / 8.0 * a_h, y_for_k3, rhs_result);
            std::array<double, RHS::N> k3 = {};
            for (int i = 0; i < RHS::N; i++) {
                k3[i] = a_h * rhs_result[i];
            }
            // k4
            std::array<double, RHS::N> y_for_k4 = {};
            for (int i = 0; i < RHS::N; i++) {
                y_for_k4[i] = a_y[i] + k4_consts[0] * k1[i] + k4_consts[1] * k2[i] + k4_consts[2] * k3[i];
            }
            m_rhs(a_t + 12.0 / 13.0 * a_h, y_for_k4, rhs_result);
            std::array<double, RHS::N> k4 = {};
            for (int i = 0; i < RHS::N; i++) {
                k4[i] = a_h * rhs_result[i];
            }
            // k5
            std::array<double, RHS::N> y_for_k5 = {};
            for (int i = 0; i < RHS::N; i++) {
                y_for_k5[i] = a_y[i] + k5_consts[0] * k1[i] + k5_consts[1] * k2[i] + k5_consts[2] * k3[i] +
                              k5_consts[3] * k4[i];
            }
            m_rhs(a_t + a_h, y_for_k5, rhs_result);
            std::array<double, RHS::N> k5 = {};
            for (int i = 0; i < RHS::N; i++) {
                k5[i] = a_h * rhs_result[i];
            }
            // k6
            std::array<double, RHS::N> y_for_k6 = {};
            for (int i = 0; i < RHS::N; i++) {
                y_for_k6[i] = a_y[i] + k6_consts[0] * k1[i] + k6_consts[1] * k2[i] + k6_consts[2] * k3[i] +
                              k6_consts[3] * k4[i] + k6_consts[4] * k5[i];
            }
            m_rhs(a_t + 1. / 2. * a_h, y_for_k6, rhs_result);
            std::array<double, RHS::N> k6 = {};
            for (int i = 0; i < RHS::N; i++) {
                k6[i] = a_h * rhs_result[i];
            }
            // delta
            std::array<double, 6> delta = {};
            for (int i = 0; i < delta.size(); i++) {
                delta[i] = delta_consts[0] * k1[i] + delta_consts[1] * k2[i] + delta_consts[2] * k3[i] +
                           delta_consts[3] * k4[i] + delta_consts[4] * k5[i] + delta_consts[5] * k6[i];
            }
            double delta_norm = 0;
            for (int i = 0; i < delta.size(); i++) {
                delta_norm += delta[i] * delta[i];
            }
            delta_norm = sqrt(delta_norm); // sqrt(delta_1^2 + ... + delta_6^2)

            if (delta_norm > eps) {
                a_h = 0.9 * a_h * pow(eps / delta_norm, 1. / 5.);
            } else {
                for (int i = 0; i < RHS::N; i++) {
                    a_y_next = a_y[i] + gamma_consts[0] * k1 + gamma_consts[1] * k2 + gamma_consts[2] * k3 +
                               gamma_consts[3] * k4 + gamma_consts[4] * k5 + gamma_consts[5] * k6;
                }
                a_t += a_h;
                flag = false;
            }
        }
    }
};

template <typename RHS>
class TimeStepper_RKF45_GSL {};

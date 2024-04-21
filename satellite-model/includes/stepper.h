#pragma once

#include <array>
#include <cstddef>
#include "function.h"

template <size_t N_steps>
class EverhartStepper {
public:
    explicit EverhartStepper(BaseFunction &rhs) : m_rhs(rhs) {
    }

    // returns <y, dy>
    std::pair<Arg, Arg> MakeStep(const double &a_h, const double &a_t, const Arg &a_y, const Arg &a_dy) {
        // Initialize approximations
        Values B_coef{}; Values div_diff{};
        Values y_values_next{}; Values dy_values_next{};
        auto [y_values, dy_values] = InitArgs(a_h, a_t, a_y, a_dy);
        for (int i = 0; i < N_steps; ++i) {
            // Get divided differences
            GetDifferences(a_h, a_t, y_values, dy_values, div_diff);
            // Convolve coefficients into B_j
            Convolve5(a_h, a_t, div_diff, B_coef);
            // Update approximation
            UpdateApproximation(
                a_h, y_values, dy_values, B_coef, y_values_next, dy_values_next
            );
            // Check for convergence
            double overall_shift = 0;
            for (int j = 0; j <= K; ++j) {
                overall_shift += (y_values_next[i] - y_values[i]).SquaredNorm();
                overall_shift += (dy_values_next[i] - dy_values[i]).SquaredNorm();
            }
            y_values = y_values_next;
            dy_values = dy_values_next;
            if (overall_shift < EPS)
                break;
        }
        return std::make_pair(y_values.back(), dy_values.back());
    }

private:
    BaseFunction &m_rhs;
    static constexpr size_t K = 5;
    static constexpr double EPS = 1e-4;

    using Values = std::array<Arg, K + 1>;

    // Evaluates delta of n-th step of uniform partition
    double EvalDelta(const double a_h, const size_t n) {
        return a_h * ((double)n / (double)K);
    }

    // Initial approximation
    std::pair<Values, Values>
    InitArgs(const double &a_h, const double &a_t, const Arg &a_y, const Arg &a_dy) {
        Values y{};
        Values dy{};
        Arg F0 = m_rhs(a_t, a_y, a_dy);
        for (int i = 0; i <= K; ++i) {
            double delta = EvalDelta(a_h, i);
            dy[i] = a_dy + F0 * delta;
            y[i] = a_y + a_dy * delta + F0 * (delta * delta / 2.0);
        }
        return std::make_pair(y, dy);
    }

    // Calculate needed divided differences
    void GetDifferences(
        double a_h,
        double a_t,
        const Values &y_values,
        const Values &dy_values,
        Values &div_diff
    ) {
        Values cur_layer{};
        for (int i = 0; i <= K; ++i) {
            double ti = a_t + EvalDelta(a_h, i);
            cur_layer[i] = m_rhs(ti, y_values[i], dy_values[i]);
        }

        Values next_layer{};
        div_diff[0] = cur_layer[0];
        for (int i = 1; i <= K; ++i) {
            for (int j = 0; i + j <= K; ++j) {
                Arg diff = cur_layer[i + j] - cur_layer[i + j - 1];
                next_layer[i + j] = diff * (1.0 / EvalDelta(a_h, i));
            }
            cur_layer = next_layer;
            div_diff[i] = cur_layer[i];
        }
    }

    // Precalculated coefficients for K = 5
    void Convolve5(const double a_h, const double a_t, const Values &div_diff, Values &conv) {
        std::array<double, K + 1> t{};
        for (int i = 0; i <= K; ++i) {
            t[i] = a_t + EvalDelta(a_h, i);
        }
        conv[0] = div_diff[0];
        conv[1] =
            div_diff[1] + div_diff[2] * (t[0] - t[1]) +
            div_diff[3] * ((t[0] - t[1]) * (t[0] - t[2])) +
            div_diff[4] * ((t[0] - t[1]) * (t[0] - t[2]) * (t[0] - t[3])) +
            div_diff[5] * ((t[0] - t[1]) * (t[0] - t[2]) * (t[0] - t[3]) * (t[0] - t[4]));
        conv[2] =
            div_diff[2] + div_diff[3] * (2.0 * t[0] - t[1] - t[2]) +
            div_diff[4] * (3.0 * t[0] * t[0] - 2.0 * t[1] * t[0] - 2.0 * t[2] * t[0] -
                           2.0 * t[3] * t[0] + t[1] * t[2] + t[1] * t[3] + t[2] * t[3]) +
            div_diff[5] * ((t[0] - t[1]) * (t[0] - t[2]) * (t[0] - t[3]) +
                           (t[0] - t[4]) * ((t[0] - t[1]) * (t[0] - t[2]) +
                                            (2.0 * t[0] - t[1] - t[2]) * (t[0] - t[3])));
        conv[3] = div_diff[3] + div_diff[4] * (3.0 * t[0] - t[1] - t[2] - t[3]) +
                  div_diff[5] *
                      (6.0 * t[0] * t[0] - 3.0 * t[1] * t[0] - 3.0 * t[2] * t[0] -
                       3.0 * t[3] * t[0] + 3.0 * t[4] * t[0] + t[1] * t[2] + t[1] * t[3] +
                       t[2] * t[3] + t[1] * t[4] + t[2] * t[4] + t[3] * t[4]);
        conv[4] = div_diff[4] + div_diff[5] * (4.0 * t[0] - t[1] - t[2] - t[3] - t[4]);
        conv[5] = div_diff[5];
    }

    // Update values' approximation
    void UpdateApproximation(
        const double a_h,
        const Values &y_values,
        const Values &dy_values,
        const Values &B,
        Values &y_values_next,
        Values &dy_values_next
    ) {
        for (int i = 1; i <= K; ++i) {
            dy_values_next[i] = dy_values[0];
            double delta = EvalDelta(a_h, i);
            double prod_delta = delta;
            y_values_next[i] = y_values[0] + dy_values[0] * delta;
            for (int j = 0; j <= K; ++j) {
                dy_values_next[i] += B[j] * (prod_delta / (double)(j + 1));
                prod_delta *= delta;
                y_values_next[i] +=
                    B[j] * (prod_delta / ((double)(j + 1) * (double)(j + 2)));
            }
        }
    }
};

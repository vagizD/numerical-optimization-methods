#pragma once

#include <array>
#include <utility>

template <typename RHS>
class TimeStepper_RKF45 {
private:
    RHS const *const m_rhs;
    constexpr static int N = RHS::N;

public:
    TimeStepper_RKF45(RHS const *a_rhs) : m_rhs(a_rhs){};

    std::pair<double, double> operator()(
        double a_t,
        double a_h,  // suggested step
        std::array<double, RHS::N> &a_y,
        std::array<double, RHS::N> &a_y_next
    ) {
        // RKF45 method
    }
};

template <typename RHS>
class TimeStepper_RKF45_GSL {};

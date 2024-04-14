#pragma once

#include <array>

template <typename Stepper, typename Observer>
class ODE_Integrator {
private:
    Stepper *const m_stepper;
    const Observer *const m_observer;

public:
    ODE_Integrator(Stepper *a_stepper, const Observer *a_observer)
        : m_stepper(a_stepper), m_observer(a_observer) {
    }

    void operator()(
        double a_t0,
        double a_tEnd,
        std::array<double, Stepper::N> &a_y0,
        std::array<double, Stepper::N> &a_yEnd
    ) {
    }
};
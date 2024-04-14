#pragma once

#include <array>
#include <string>

template <typename Stepper, typename Observer>
class ODE_Integrator {
private:
    Stepper *m_stepper;
    Observer *m_observer;
    std::string &m_path;
    bool m_verbose;
    bool m_save;

public:
    ODE_Integrator(
        Stepper *a_stepper,
        Observer *a_observer,
        std::string &a_path,
        bool a_verbose,
        bool a_save
    )
        : m_stepper(a_stepper),
          m_observer(a_observer),
          m_path(a_path),
          m_verbose(a_verbose),
          m_save(a_save) {
    }

    void operator()(
        double a_t0,
        double a_h,
        std::array<double, Stepper::N> &a_y0,
        std::array<double, Stepper::N> &a_yEnd
    ) {
        double t = a_t0;
        double step = a_h;
        std::array<double, Stepper::N> y = a_y0;
        std::array<double, Stepper::N> y_next = a_yEnd;

        while (true) {
            m_stepper->make_step(t, step, y, y_next);
            if (!m_observer->make_decision(t, y, m_verbose)) {
                if (m_save) {
                    m_observer->save_projectile(m_path);
                }
                a_yEnd = y;
                break;
            }
            y = y_next;
        }
    }
};

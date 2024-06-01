#pragma once
#include "function.h"

template <typename Stepper, typename Observer>
class ODE_Integrator2 {
private:
    Stepper *m_stepper;
    Observer *m_observer;
    std::string &m_path;
    bool m_verbose;
    bool m_save;

public:
    ODE_Integrator2(
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

    std::pair<Arg, Arg> operator()(double a_t0, double a_h, Arg &a_y0, Arg &a_dy0) {
        double t = a_t0;
        double step = a_h;
        Arg y = a_y0;
        Arg dy = a_dy0;

        while (true) {
            auto [y_next, dy_next] = m_stepper->MakeStep(step, t, y, dy);
            if (!m_observer->make_decision(t, y, dy, m_verbose)) {
                if (m_save) {
                    m_observer->save_projectile(m_path);
                }
                y_next = y;
                dy_next = dy;
                return std::make_pair(y_next, dy_next);
            }
            y = y_next;
            dy = dy_next;
        }
        // return std::make_pair(y_next, dy_next);
    }
};
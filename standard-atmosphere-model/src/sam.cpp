#include "sam.h"
#include <cassert>
#include <cmath>
#include "constants.h"

int getLayer(double height) {
    assert(0 <= height && height <= h.back());
    if (height <= h[1]) {
        return 0;
    } else if (height <= h[2]) {
        return 1;
    } else if (height <= h[3]) {
        return 2;
    }
    return 3;
}

Pressure::Pressure(double p0, double t0) {  // initial pressure and absolute temperature
    t[0] = t0, p[0] = p0;
    for (int i = 1; i < 4; ++i) {
        t[i] = t[i - 1] - r[i - 1] * (h[i] - h[i - 1]);
        p[i] = p[i - 1] * std::exp(
                              g / (R * r[i - 1]) *
                              std::log(1 - r[i - 1] * (h[i] - h[i - 1]) / t[i - 1])
                          );
    }
}

double Pressure::operator()(double height) const {
    int l = getLayer(height);
    if (r[l] != 0) {
        return p[l] *
               std::exp(g / (R * r[l]) * std::log(1 - r[l] * (height - h[l]) / t[l]));
    }
    return p[l] * std::exp(g * (h[l] - height) / (R * t[l]));
}

double Density::operator()(double height) const {
    int l = getLayer(height);
    return m_pressure(height) / (R * (m_pressure.t[l] - r[l] * (height - h[l])));
}

double AirSpeed::operator()(double height) const {
    return std::sqrt(m_density.m_pressure(height) / m_density(height));
}

double DragCoef::operator()(double mach) const {
    assert(mach >= 0);
    if (mach >= 2.5) {
        return m_coefs.back();
    }

    int i = static_cast<int>(mach / m_step);
    return m_coefs[i] + ((m_coefs[i + 1] - m_coefs[i]) / m_step) * (mach - i * m_step);
}

double AerodynamicDragForce::operator()(double mach, double y, double velocity, double s)
    const {
    return m_cd(mach) * m_density(y) * s * velocity * velocity / 2;
}

void RHS_Projectile::operator()(
    const double a_t,  // unused, autonomous ODE
    std::array<double, 4> &a_y,
    std::array<double, 4> &a_y_next
) const {
    double s = M_PI * m_diameter * m_diameter / 4;
    double velocity = std::sqrt(a_y[1] * a_y[1] + a_y[3] * a_y[3]);
    double mach = m_airSpeed(a_y[2]);
    double qdf = m_q(mach, a_y[2], velocity, s);

    a_y_next[0] = a_y[1];
    a_y_next[1] = -(qdf * a_y[1]) / (m_mass * velocity);
    a_y_next[2] = a_y[3];
    a_y_next[3] = -(qdf * a_y[3]) / (m_mass * velocity) - g;
}

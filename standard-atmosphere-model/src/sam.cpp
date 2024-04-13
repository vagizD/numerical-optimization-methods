#include "sam.h"
#include <cassert>
#include <cmath>
#include "constants.h"

Pressure::Pressure(double p0, double t0) {
    t[0] = t0, p[0] = p0;
    for (int i = 1; i < 4; ++i) {
        t[i] = t[i - 1] - r[i - 1] * (h[i] - h[i - 1]);
        p[i] = p[i - 1] * std::exp(
                              g / (R * r[i - 1]) *
                              std::log(1 - r[i - 1] * (h[i] - h[i - 1]) / t[i - 1])
                          );
    }
}

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

double Pressure::operator()(double height) {
    int l = getLayer(height);
    if (r[l] != 0) {
        return p[l] *
               std::exp(g / (R * r[l]) * std::log(1 - r[l] * (height - h[l]) / t[l]));
    }
    return p[l] * std::exp(g * (h[l] - height) / (R * t[l]));
}

double Density::operator()(double height) {
    int l = getLayer(height);
    return pressure(height) / (R * (pressure.t[l] - r[l] * (height - h[l])));
}

double AirSpeed::operator()(double height) {
    return std::sqrt(density.pressure(height) / density(height));
}

double DragCoef::operator()(double mach) {
    assert(mach >= 0);
    if (mach >= 2.5) {
        return coefs.back();
    }

    int i = static_cast<int>(mach / step);
    return coefs[i] + ((coefs[i + 1] - coefs[i]) / step) * (mach - i * step);
}

double
AerodynamicDragForce::operator()(double mach, double y, double velocity, double s) {
    return cd(mach) * density(y) * s * velocity * velocity / 2;
}

void nextState(
    AerodynamicDragForce q,
    double mass,
    double mach,
    double s,
    const std::array<double, 4> &u,
    std::array<double, 4> &result
) {
    double velocity = std::sqrt(u[1] * u[1] + u[3] * u[3]);
    double qdf = q(mach, u[2], velocity, s);
    result[0] = u[1];
    result[1] = -(qdf * u[1]) / (mass * velocity);
    result[2] = u[3];
    result[3] = -(qdf * u[3]) / (mass * velocity) - g;
}

void RHS_Velocity::operator()(
    const double a_t,
    std::array<double, 4> &a_v,
    std::array<double, 4> &result
) {
    Pressure pressure(m_p0, m_t0);
    Density density(pressure);
    AirSpeed airSpeed(density);
    DragCoef cd;
    AerodynamicDragForce Q(cd, density);
    nextState(Q, m_mass, airSpeed(a_v[2]), m_s, a_v, result);
}

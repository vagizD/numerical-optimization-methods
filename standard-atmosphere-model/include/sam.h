#pragma once

#include <array>
#include <cmath>

class Density;
class AirSpeed;

class Pressure {
public:
    Pressure(double p0, double t0);
    double operator()(double height) const;

private:
    friend class Density;
    std::array<double, 4> p;
    std::array<double, 5> t;
};

class Density {
public:
    explicit Density(Pressure p) : m_pressure(p){};
    double operator()(double height) const;

private:
    friend class AirSpeed;
    Pressure m_pressure;
};

class AirSpeed {
public:
    explicit AirSpeed(Density d) : m_density(d){};
    double operator()(double height) const;

private:
    Density m_density;
};

class DragCoef {
public:
    DragCoef(){};
    double operator()(double mach) const;

private:
    double m_step = 0.1;
    std::array<double, 26> m_coefs = {
        0.00,  // 0.0
        0.10,  // 0.1
        0.10,  // 0.2
        0.10,  // 0.3
        0.10,  // 0.4
        0.10,  // 0.5
        0.10,  // 0.6
        0.11,  // 0.7
        0.12,  // 0.8
        0.17,  // 0.9
        0.38,  // 1.0
        0.34,  // 1.1
        0.33,  // 1.2
        0.32,  // 1.3
        0.31,  // 1.4
        0.30,  // 1.5
        0.30,  // 1.6
        0.29,  // 1.7
        0.27,  // 1.8
        0.26,  // 1.9
        0.25,  // 2.0
        0.24,  // 2.1
        0.23,  // 2.2
        0.22,  // 2.3
        0.21,  // 2.4
        0.20,  // 2.5
    };
};

class AerodynamicDragForce {
public:
    AerodynamicDragForce(DragCoef dragCoef, Density d) : m_cd(dragCoef), m_density(d){};
    double operator()(double mach, double y, double velocity, double s) const;

private:
    DragCoef m_cd;
    Density m_density;
};

class RHS_Projectile {
private:
    double m_mass;
    double m_diameter;

    const Pressure m_pressure;
    const Density m_density;
    const AirSpeed m_airSpeed;
    DragCoef m_cd;
    AerodynamicDragForce m_q;

public:
    constexpr static int N = 4;
    RHS_Projectile(
        double a_pressure0,
        double a_temperature0,
        double a_mass,
        double a_diameter
    )
        : m_mass(a_mass),
          m_diameter(a_diameter),
          m_pressure(a_pressure0, a_temperature0),
          m_density(m_pressure),
          m_airSpeed(m_density),
          m_cd(),
          m_q(m_cd, m_density){};
    void operator()(
        double a_t,
        const std::array<double, N> &a_y,
        std::array<double, N> &a_y_next
    ) const;
};

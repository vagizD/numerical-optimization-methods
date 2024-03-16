#pragma once

#include <array>

class Density;
class AirSpeed;

class Pressure {
public:
    Pressure(double p0, double t0);
    double operator()(double height);
private:
    friend class Density;
    std::array<double, 4> p;
    std::array<double, 5> t;
};

class Density {
    Density(Pressure p): pressure(p) {};
    double operator()(double height);
private:
    friend class AirSpeed;
    Pressure pressure;
};

class AirSpeed {
    AirSpeed(Density d): density(d) {};
    double operator()(double height);
private:
    Density density;
};

class DragCoef {
    DragCoef() {};
    double operator()(double mach);
private:
    double step = 0.1;
    std::array<double, 26> coefs = {
            0.00,      // 0.0
            0.10,      // 0.1
            0.10,      // 0.2
            0.10,      // 0.3
            0.10,      // 0.4
            0.10,      // 0.5
            0.10,      // 0.6
            0.11,      // 0.7
            0.12,      // 0.8
            0.17,      // 0.9
            0.38,      // 1.0
            0.34,      // 1.1
            0.33,      // 1.2
            0.32,      // 1.3
            0.31,      // 1.4
            0.30,      // 1.5
            0.30,      // 1.6
            0.29,      // 1.7
            0.27,      // 1.8
            0.26,      // 1.9
            0.25,      // 2.0
            0.24,      // 2.1
            0.23,      // 2.2
            0.22,      // 2.3
            0.21,      // 2.4
            0.20,      // 2.5
    };
};
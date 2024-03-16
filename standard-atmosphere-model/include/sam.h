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

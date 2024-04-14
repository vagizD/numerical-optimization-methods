#include <cmath>
#include <iostream>
#include "constants.h"
#include "sam.h"

int main() {
    double p0 = 101325;
    double t0 = 288.15;
    double mass = 106.0;
    double diameter = 0.216;

    double alpha = 30;  // degrees
    double v0 = 1640;   // meter/second

    std::array<double, 4> u = {
        10.0, v0 * std::cos(alpha * M_PI / 180), 5.0, v0 * std::sin(alpha * M_PI / 180)
    };
    std::array<double, 4> u_diff = {};

    RHS_Velocity rhs = RHS_Velocity(p0, t0, mass, diameter);
    rhs(t0, u, u_diff);

    std::cout << "u : "
              << "< " << u[0] << ", " << u[1] << ", " << u[2] << ", " << u[3] << " >"
              << std::endl;
    std::cout << "u': "
              << "< " << u_diff[0] << ", " << u_diff[1] << ", " << u_diff[2] << ", "
              << u_diff[3] << " >" << std::endl;

    return 0;
}
#include <cmath>
#include <iostream>
#include "constants.h"
#include "sam.h"

int main() {
    double p0 = 101325;
    double t0 = 288.15;
    double mass = 106.0;
    double diameter = 0.216;
    double s = M_PI * diameter * diameter / 4;

    std::array<double, 4> u = {
        10.0, 1483.104, 5.0, 900.05};  // sqrt(v_x^2 + v_y^2) ~ 1640.0
    std::array<double, 4> u_diff = {};

    Pressure pressure(p0, t0);
    Density density(pressure);
    AirSpeed airSpeed(density);
    DragCoef cd;
    AerodynamicDragForce Q(cd, density);
    nextState(Q, mass, airSpeed(u[0]), s, u, u_diff);

    std::cout << "u : "
              << "< " << u[0] << " " << u[1] << " " << u[2] << " " << u[3] << " >"
              << std::endl;
    std::cout << "u': "
              << "< " << u_diff[0] << " " << u_diff[1] << " " << u_diff[2] << " "
              << u_diff[3] << " >" << std::endl;

    return 0;
}
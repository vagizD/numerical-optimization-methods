#include <cmath>
#include <string>
#include "constants.h"
#include "integrator.h"
#include "observer.h"
#include "sam.h"
#include "timestepper.h"

int main() {
    double p0 = 101325;
    double t0 = 288.15;
    double mass = 106.0;
    double diameter = 0.216;

    double alpha = 5;  // degrees
    double v0 = 1640;  // meter/second
    double x0 = 10;
    double y0 = 0;

    std::string path = "../projectiles/trajectory/trajectory5.csv";

    auto rhs = RHS_Projectile(p0, t0, mass, diameter);
    auto stepper = TimeStepper_RKF45<RHS_Projectile>(&rhs);
    auto observer = SimpleObserver(x0, y0);
    auto integrator = ODE_Integrator<TimeStepper_RKF45<RHS_Projectile>, SimpleObserver>(
        &stepper, &observer, path, false, false
    );

    double t_start = 0;
    double h_start = 1e-3;
    std::array<double, 4> u = {
        x0, v0 * std::cos(alpha * M_PI / 180), y0, v0 * std::sin(alpha * M_PI / 180)
    };
    std::array<double, 4> result = {};

    integrator(t_start, h_start, u, result);

    return 0;
}

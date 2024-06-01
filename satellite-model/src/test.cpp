#include <cmath>
#include <string>
#include "integrator.h"
#include "observer.h"
#include "stepper.h"

struct TestFunction : BaseFunction {
private:
    const double mu = 398600.4;
    const double J = 0.0010827;
    const double R = 6378.137;

public:
    Arg operator()(const double a_t, const Arg &a_y, const Arg &a_dy) override {
        Arg new_arg;
        double x = a_y[0];
        double y = a_y[1];
        double z = a_y[2];
        new_arg[0] = -mu * a_y[0] / pow(a_y.SquaredNorm(), 3) -
                     J * R * R * mu *
                         (-7.5 * a_y[2] * a_y[2] * a_y[0] / pow(a_y.SquaredNorm(), 7) +
                          1.5 * a_y[0] / pow(a_y.SquaredNorm(), 5));
        new_arg[1] = -mu * a_y[1] / pow(a_y.SquaredNorm(), 3) -
                     J * R * R * mu *
                         (-7.5 * a_y[2] * a_y[2] * a_y[1] / pow(a_y.SquaredNorm(), 7) +
                          1.5 * a_y[1] / pow(a_y.SquaredNorm(), 5));
        new_arg[2] =
            -mu * a_y[2] / pow(a_y.SquaredNorm(), 3) -
            J * R * R * mu *
                (1.5 * a_y[2] *
                     (2 * a_y[0] * a_y[0] + 2 * a_y[1] * a_y[1] - 3 * a_y[2] * a_y[2]) /
                     pow(a_y.SquaredNorm(), 7) +
                 1.5 * a_y[2] / pow(a_y.SquaredNorm(), 5));
        return new_arg;
    }
};

int main() {
    double t_start = 0.0;
    double t_end = 3.1e7;
    double h_start = 1.0;
    double orbit = 45000;
    double R = 6378.137;

    TestFunction func;
    EverhartStepper<50> stepper(func);

    double a = 7500;
    double mu = 398600.4;
    double v0 = sqrt(mu / a);
    Arg test_y = Arg(0, 0, a);
    Arg test_dy = Arg(v0, 0, 0);
    auto observer = SatelliteObserver(0, 0, a, orbit, R, t_end);
    std::string path = "../satellite/trajectory/trajectory_1.csv";
    auto integrator = ODE_Integrator2<EverhartStepper<50>, SatelliteObserver>(
        &stepper, &observer, path, true, true
    );

    auto [y_res, dy_res] = integrator(t_start, h_start, test_y, test_dy);

    return 0;
}

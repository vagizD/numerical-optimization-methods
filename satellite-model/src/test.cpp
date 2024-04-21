#include "stepper.h"

struct TestFunction : BaseFunction {
    Arg operator()(const double a_t, const Arg &a_y, const Arg &a_dy) override {
        return a_y * a_t, a_y * a_dy, a_y;
    }
};

int main() {
    TestFunction func;
    EverhartStepper<10> stepper(func);

    double test_t = 0;
    double test_h = 1e-3;
    Arg test_y = Arg(1, 2, 3);
    Arg test_dy = Arg(3, 2, 1);

    auto [next_y, next_dy] = stepper.MakeStep(test_h, test_t, test_y, test_dy);
    std::cout << next_y << next_dy;

    //    3.23993e-14 6.47653e-14 9.71314e-14
    //    1.37807e-10 2.75447e-10 4.13087e-10

    return 0;
}
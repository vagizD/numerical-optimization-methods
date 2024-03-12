#include "tests.h"
#include <functional>
#include <iostream>
#include "aad.h"
#include "differentiator.h"

double F(double x, double y) {
    return std::cos(x * 5) / (x * x + y * y);
}

AAD22 autoF(AAD22 x, AAD22 y) {
    return cos(x * 5) / (x * x + y * y);
}

double dFy(double x, double y) {
    return -(y * 2 * std::cos(x * 5)) / ((x * x + y * y) * (x * x + y * y));
}

int main() {
    {
        double l_x = -50, r_x = 50, step_x = 0.25;
        double l_y = 1, r_y = 100, step_y = 0.25;
        double err_s3 = makeTests<Derivative::Y, DiffMethod::Stencil3, double>(
            F, dFy, l_x, r_x, step_x, l_y, r_y, step_y
        );
        double err_s3e = makeTests<Derivative::Y, DiffMethod::Stencil3Extra, double>(
            F, dFy, l_x, r_x, step_x, l_y, r_y, step_y
        );
        double err_s5 = makeTests<Derivative::Y, DiffMethod::Stencil5, double>(
            F, dFy, l_x, r_x, step_x, l_y, r_y, step_y
        );
        double err_s5e = makeTests<Derivative::Y, DiffMethod::Stencil5Extra, double>(
            F, dFy, l_x, r_x, step_x, l_y, r_y, step_y
        );
        double err_auto = makeTests<Derivative::Y, DiffMethod::FwdADD, AAD22>(
            autoF, dFy, l_x, r_x, step_x, l_y, r_y, step_y
        );
        std::cout
            << "... TESTING F = cos(5x) / (x^2 + y^2), (x, y) ∈ [-50, 50] x [1, 100]"
            << std::endl;
        std::cout << "=>  STENCIL3     : " << err_s3 << std::endl;
        std::cout << "=>  STENCIL3EXTRA: " << err_s3e << std::endl;
        std::cout << "=>  STENCIL5     : " << err_s5 << std::endl;
        std::cout << "=>  STENCIL5EXTRA: " << err_s5e << std::endl;
        std::cout << "=>  AAD          : " << err_auto << std::endl;
        std::cout << std::endl;
    }
    {
        auto f = [](double x, double y) { return std::exp(std::sin(x * y) + 1) * 3; };

        auto af = [](AAD22 x, AAD22 y) { return exp(sin(x * y) + 1) * 3; };

        auto dfx = [](double x, double y) {
            return y * 3 * std::exp(std::sin(x * y) + 1) * std::cos(x * y);
        };

        double l_x = -10, r_x = 10, step_x = 0.1;
        double l_y = -10, r_y = 10, step_y = 0.1;
        double err_s3 = makeTests<Derivative::X, DiffMethod::Stencil3, double>(
            f, dfx, l_x, r_x, step_x, l_y, r_y, step_y
        );
        double err_s3e = makeTests<Derivative::X, DiffMethod::Stencil3Extra, double>(
            f, dfx, l_x, r_x, step_x, l_y, r_y, step_y
        );
        double err_s5 = makeTests<Derivative::X, DiffMethod::Stencil5, double>(
            f, dfx, l_x, r_x, step_x, l_y, r_y, step_y
        );
        double err_s5e = makeTests<Derivative::X, DiffMethod::Stencil5Extra, double>(
            f, dfx, l_x, r_x, step_x, l_y, r_y, step_y
        );
        double err_auto = makeTests<Derivative::X, DiffMethod::FwdADD, AAD22>(
            af, dfx, l_x, r_x, step_x, l_y, r_y, step_y
        );
        std::cout
            << "... TESTING F = 3 * exp(sin(xy) + 1), (x, y) ∈ [-10, 10] x [-10, 10]"
            << std::endl;
        std::cout << "=>  STENCIL3     : " << err_s3 << std::endl;
        std::cout << "=>  STENCIL3EXTRA: " << err_s3e << std::endl;
        std::cout << "=>  STENCIL5     : " << err_s5 << std::endl;
        std::cout << "=>  STENCIL5EXTRA: " << err_s5e << std::endl;
        std::cout << "=>  AAD          : " << err_auto << std::endl;
        std::cout << std::endl;
    }
    {
        struct f {
            double operator()(double x, double y) {
                return std::sin(x * y) / std::exp(x - y + 1);
            }
        } f;

        struct af {
            AAD22 operator()(AAD22 x, AAD22 y) {
                return sin(x * y) / exp(x - y + 1);
            }
        } af;

        struct dfxx {
            double operator()(double x, double y) {
                return -std::exp(y - x - 1) *
                       ((y * y - 1) * std::sin(x * y) + 2 * y * std::cos(x * y));
            }
        } dfxx;

        double l_x = -2, r_x = 2, step_x = 0.1;
        double l_y = -2, r_y = 2, step_y = 0.1;
        double err_s3 = makeTests<Derivative::XX, DiffMethod::Stencil3, double>(
            f, dfxx, l_x, r_x, step_x, l_y, r_y, step_y
        );
        double err_s3e = makeTests<Derivative::XX, DiffMethod::Stencil3Extra, double>(
            f, dfxx, l_x, r_x, step_x, l_y, r_y, step_y
        );
        double err_s5 = makeTests<Derivative::XX, DiffMethod::Stencil5, double>(
            f, dfxx, l_x, r_x, step_x, l_y, r_y, step_y
        );
        double err_s5e = makeTests<Derivative::XX, DiffMethod::Stencil5Extra, double>(
            f, dfxx, l_x, r_x, step_x, l_y, r_y, step_y
        );
        double err_auto = makeTests<Derivative::XX, DiffMethod::FwdADD, AAD22>(
            af, dfxx, l_x, r_x, step_x, l_y, r_y, step_y
        );
        std::cout
            << "... TESTING F = sin(xy) / exp(x - y + 1), (x, y) ∈ [-2, 2] x [-2, 2]"
            << std::endl;
        std::cout << "=>  STENCIL3     : " << err_s3 << std::endl;
        std::cout << "=>  STENCIL3EXTRA: " << err_s3e << std::endl;
        std::cout << "=>  STENCIL5     : " << err_s5 << std::endl;
        std::cout << "=>  STENCIL5EXTRA: " << err_s5e << std::endl;
        std::cout << "=>  AAD          : " << err_auto << std::endl;
        std::cout << std::endl;
    }
    {
        auto f = [](double x, double y) {
            return std::sin(x + y + M_PI) * std::cos(x - y);
        };

        auto af = [](AAD22 x, AAD22 y) { return sin(x + y + M_PI) * cos(x - y); };

        auto dfyy = [](double x, double y) {
            return 2 * (std::sin(x + y) * std::cos(x - y) -
                        std::sin(x - y) * std::cos(x + y));
        };

        double l_x = -10, r_x = 10, step_x = 0.1;
        double l_y = -10, r_y = 10, step_y = 0.1;
        double err_s3 = makeTests<Derivative::YY, DiffMethod::Stencil3, double>(
            f, dfyy, l_x, r_x, step_x, l_y, r_y, step_y
        );
        double err_s3e = makeTests<Derivative::YY, DiffMethod::Stencil3Extra, double>(
            f, dfyy, l_x, r_x, step_x, l_y, r_y, step_y
        );
        double err_s5 = makeTests<Derivative::YY, DiffMethod::Stencil5, double>(
            f, dfyy, l_x, r_x, step_x, l_y, r_y, step_y
        );
        double err_s5e = makeTests<Derivative::YY, DiffMethod::Stencil5Extra, double>(
            f, dfyy, l_x, r_x, step_x, l_y, r_y, step_y
        );
        double err_auto = makeTests<Derivative::YY, DiffMethod::FwdADD, AAD22>(
            af, dfyy, l_x, r_x, step_x, l_y, r_y, step_y
        );
        std::cout << "... TESTING F = sin(x + y + π) * cos(x - y), (x, y) ∈ [-10, 10] x "
                     "[-10, 10]"
                  << std::endl;
        std::cout << "=>  STENCIL3     : " << err_s3 << std::endl;
        std::cout << "=>  STENCIL3EXTRA: " << err_s3e << std::endl;
        std::cout << "=>  STENCIL5     : " << err_s5 << std::endl;
        std::cout << "=>  STENCIL5EXTRA: " << err_s5e << std::endl;
        std::cout << "=>  AAD          : " << err_auto << std::endl;
        std::cout << std::endl;
    }
    {
        auto f = [](double x, double y) { return std::exp(x / y) * std::sin(x * 5); };

        auto af = [](AAD22 x, AAD22 y) { return exp(x / y) * sin(x * 5); };

        auto dfxy = [](double x, double y) {
            return -(std::exp(x / y) * (x * std::sin(x * 5) +
                                        y * (std::sin(x * 5) + 5 * x * std::cos(x * 5)))
                   ) /
                   (y * y * y);
        };

        double l_x = -5, r_x = 5, step_x = 0.1;
        double l_y = 1, r_y = 5, step_y = 0.1;
        double err_s3 = makeTests<Derivative::XY, DiffMethod::Stencil3, double>(
            f, dfxy, l_x, r_x, step_x, l_y, r_y, step_y
        );
        double err_s3e = makeTests<Derivative::XY, DiffMethod::Stencil3Extra, double>(
            f, dfxy, l_x, r_x, step_x, l_y, r_y, step_y
        );
        double err_s5 = makeTests<Derivative::XY, DiffMethod::Stencil5, double>(
            f, dfxy, l_x, r_x, step_x, l_y, r_y, step_y
        );
        double err_s5e = makeTests<Derivative::XY, DiffMethod::Stencil5Extra, double>(
            f, dfxy, l_x, r_x, step_x, l_y, r_y, step_y
        );
        double err_auto = makeTests<Derivative::XY, DiffMethod::FwdADD, AAD22>(
            af, dfxy, l_x, r_x, step_x, l_y, r_y, step_y
        );
        std::cout << "... TESTING F = exp(x / y) * sin(5x), (x, y) ∈ [-5, 5] x [1, 5]"
                  << std::endl;
        std::cout << "=>  STENCIL3     : " << err_s3 << std::endl;
        std::cout << "=>  STENCIL3EXTRA: " << err_s3e << std::endl;
        std::cout << "=>  STENCIL5     : " << err_s5 << std::endl;
        std::cout << "=>  STENCIL5EXTRA: " << err_s5e << std::endl;
        std::cout << "=>  AAD          : " << err_auto << std::endl;
        std::cout << std::endl;
    }

    // LOCAL RESULTS
    // ... TESTING F = cos(5x) / (x^2 + y^2), (x, y) ∈ [-50, 50] x [1, 100]
    // =>  STENCIL3     : 3.99989e-08
    // =>  STENCIL3EXTRA: 3.82239e-12
    // =>  STENCIL5     : 1.51568e-12
    // =>  STENCIL5EXTRA: 5.78493e-12
    // =>  AAD          : 0
    //
    // ... TESTING F = 3 * exp(sin(xy) + 1), (x, y) ∈ [-10, 10] x [-10, 10]
    // =>  STENCIL3     : 0.00464276
    // =>  STENCIL3EXTRA: 2.16268e-09
    // =>  STENCIL5     : 4.89897e-07
    // =>  STENCIL5EXTRA: 5.99881e-09
    // =>  AAD          : 4.26326e-14
    //
    // ... TESTING F = sin(xy) / exp(x - y + 1), (x, y) ∈ [-2, 2] x [-2, 2]
    // =>  STENCIL3     : 1.14418e-06
    // =>  STENCIL3EXTRA: 4.10549e-05
    // =>  STENCIL5     : 6.04692e-07
    // =>  STENCIL5EXTRA: 5.30169e-05
    // =>  AAD          : 2.13163e-14
    //
    // ... TESTING F = sin(x + y + π) * cos(x - y), (x, y) ∈ [-10, 10] x [-10, 10]
    // =>  STENCIL3     : 6.11949e-07
    // =>  STENCIL3EXTRA: 1.96232e-05
    // =>  STENCIL5     : 2.48412e-07
    // =>  STENCIL5EXTRA: 2.45291e-05
    // =>  AAD          : 4.10783e-15
    //
    // ... TESTING F = exp(x / y) * sin(5x), (x, y) ∈ [-5, 5] x [1, 5]
    // =>  STENCIL3     : 0.00246631
    // =>  STENCIL3EXTRA: 4.29148e-05
    // =>  STENCIL5     : 5.11332e-07
    // =>  STENCIL5EXTRA: 8.88344e-05
    // =>  AAD          : 4.54747e-13

    return 0;
}

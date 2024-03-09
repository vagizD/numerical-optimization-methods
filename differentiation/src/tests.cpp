#include <functional>
#include <iostream>
#include "aad.h"
#include "differentiator.h"

double f1TestStencil(double x, double y) {
    return x * y / (x + y);
}

AAD22 f1TestAAD(AAD22 x, AAD22 y) {
    return x * y / (x + y);
}

int main() {
    {
        double da = Differentiator<
            Derivative::X, DiffMethod::FwdADD, std::function<AAD22(AAD22, AAD22)>>(
            f1TestAAD, 1, 1
        );
        double ds3 = Differentiator<
            Derivative::X, DiffMethod::Stencil3, std::function<double(double, double)>>(
            f1TestStencil, 1, 1
        );
        double ds3e = Differentiator<
            Derivative::X, DiffMethod::Stencil3Extra,
            std::function<double(double, double)>>(f1TestStencil, 1, 1);
        double ds5 = Differentiator<
            Derivative::X, DiffMethod::Stencil5, std::function<double(double, double)>>(
            f1TestStencil, 1, 1
        );
        double ds5e = Differentiator<
            Derivative::X, DiffMethod::Stencil5Extra,
            std::function<double(double, double)>>(f1TestStencil, 1, 1);
        std::cout << "Testing for F(x, y) = xy / (x + y)" << std::endl;
        std::cout << "AAD:           F'x(1, 1) = " << da << std::endl;
        std::cout << "Stencil3:      F'x(1, 1) = " << ds3 << std::endl;
        std::cout << "Stencil3Extra: F'x(1, 1) = " << ds3e << std::endl;
        std::cout << "Stencil5:      F'x(1, 1) = " << ds5 << std::endl;
        std::cout << "Stencil5Extra: F'x(1, 1) = " << ds5e << std::endl;
        std::cout << std::endl;
    }
    {
        auto f = [](double x, double y) { return -(std::sin(-x) * std::exp(y * 2)); };
        auto g = [](AAD22 x, AAD22 y) { return -(sin(-x) * exp(y * 2)); };

        double da = Differentiator<
            Derivative::XY, DiffMethod::FwdADD, std::function<AAD22(AAD22, AAD22)>>(
            g, M_PI, 0
        );
        double ds3 = Differentiator<
            Derivative::XY, DiffMethod::Stencil3, std::function<double(double, double)>>(
            f, M_PI, 0
        );
        double ds3e = Differentiator<
            Derivative::XY, DiffMethod::Stencil3Extra,
            std::function<double(double, double)>>(f, M_PI, 0);
        double ds5 = Differentiator<
            Derivative::XY, DiffMethod::Stencil5, std::function<double(double, double)>>(
            f, M_PI, 0
        );
        double ds5e = Differentiator<
            Derivative::XY, DiffMethod::Stencil5Extra,
            std::function<double(double, double)>>(f, M_PI, 0);
        std::cout << "Testing for F(x, y) = -(sin(-x) * exp(y))" << std::endl;
        std::cout << "AAD:           F'xy(π, 0) = " << da << std::endl;
        std::cout << "Stencil3:      F'xy(π, 0) = " << ds3 << std::endl;
        std::cout << "Stencil3Extra: F'xy(π, 0) = " << ds3e << std::endl;
        std::cout << "Stencil5:      F'xy(π, 0) = " << ds5 << std::endl;
        std::cout << "Stencil5Extra: F'xy(π, 0) = " << ds5e << std::endl;
        std::cout << std::endl;
    }
    {
        struct F {
            double operator()(double x, double y) {
                return std::sin(std::exp(x + y)) / std::cos(std::exp(x + y));
            }
        } F;

        struct G {
            AAD22 operator()(AAD22 x, AAD22 y) {
                return sin(exp(x + y)) / cos(exp(x + y));
            }
        } G;

        double X = std::log(std::sqrt(M_PI));
        double Y = std::log(std::sqrt(M_PI));
        double da = Differentiator<
            Derivative::Y, DiffMethod::FwdADD, std::function<AAD22(AAD22, AAD22)>>(
            G, X, Y
        );
        double ds3 = Differentiator<
            Derivative::Y, DiffMethod::Stencil3, std::function<double(double, double)>>(
            F, X, Y
        );
        double ds3e = Differentiator<
            Derivative::Y, DiffMethod::Stencil3Extra,
            std::function<double(double, double)>>(F, X, Y);
        double ds5 = Differentiator<
            Derivative::Y, DiffMethod::Stencil5, std::function<double(double, double)>>(
            F, X, Y
        );
        double ds5e = Differentiator<
            Derivative::Y, DiffMethod::Stencil5Extra,
            std::function<double(double, double)>>(F, X, Y);
        std::cout << "Testing for F(x, y) = tg(exp(x + y))" << std::endl;
        std::cout << "AAD:           F'y(ln(π)/2, ln(π)/2) = " << da << std::endl;
        std::cout << "Stencil3:      F'y(ln(π)/2, ln(π)/2) = " << ds3 << std::endl;
        std::cout << "Stencil3Extra: F'y(ln(π)/2, ln(π)/2) = " << ds3e << std::endl;
        std::cout << "Stencil5:      F'y(ln(π)/2, ln(π)/2) = " << ds5 << std::endl;
        std::cout << "Stencil5Extra: F'y(ln(π)/2, ln(π)/2) = " << ds5e << std::endl;
        std::cout << std::endl;
    }
    {
        auto f = [](double x, double y) {
            return std::sin(y * y + M_PI * std::cos(-std::exp((x + 1) / (y * 3)))) *
                       std::exp(std::cos(x - y) / std::sin(x - y)) +
                   1;
        };

        AAD22 X = AAD22(Variable::X, 1);
        AAD22 Y = AAD22(Variable::Y, 2);
        AAD22 A = Y * Y;
        A += cos(-exp((X + 1) / (Y * 3))) * M_PI;
        A = sin(A);
        AAD22 B = cos(X - Y);
        B /= sin(X - Y);
        A *= exp(B);
        A += 1;

        double da = A.get_derivative(Derivative::XX);
        double ds3 = Differentiator<
            Derivative::XX, DiffMethod::Stencil3, std::function<double(double, double)>>(
            f, 1, 2
        );
        double ds3e = Differentiator<
            Derivative::XX, DiffMethod::Stencil3Extra,
            std::function<double(double, double)>>(f, 1, 2);
        double ds5 = Differentiator<
            Derivative::XX, DiffMethod::Stencil5, std::function<double(double, double)>>(
            f, 1, 2
        );
        double ds5e = Differentiator<
            Derivative::XX, DiffMethod::Stencil5Extra,
            std::function<double(double, double)>>(f, 1, 2);
        std::cout << "Testing for F(x, y) = :)" << std::endl;
        std::cout << "AAD:           F'xx(1, 2) = " << da << std::endl;
        std::cout << "Stencil3:      F'xx(1, 2) = " << ds3 << std::endl;
        std::cout << "Stencil3Extra: F'xx(1, 2) = " << ds3e << std::endl;
        std::cout << "Stencil5:      F'xx(1, 2) = " << ds5 << std::endl;
        std::cout << "Stencil5Extra: F'xx(1, 2) = " << ds5e << std::endl;
        std::cout << std::endl;
    }

    return 0;
}
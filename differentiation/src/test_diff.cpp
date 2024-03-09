#include "../include/diff.h"

#include <cmath>
#include <iostream>

int main() {
    AAD22 X = AAD22(Variable::X, M_PI);
    AAD22 Y = AAD22(Variable::Y, M_PI);
    AAD22 Z = AAD22(Variable::Y, M_PI);
    AAD22 F = sin(X);
    AAD22 G = cos(Y);
    AAD22 H = exp(Z);

    std::cout << "F(x, y) = sin(x)" << std::endl;
    std::cout << "F(π, *) = " << F.get_value() << std::endl;
    std::cout << "F'x(π, *) = " << F.get_derivative(Derivative::X) << std::endl;
    std::cout << "F'y(π, *) = " << F.get_derivative(Derivative::Y) << std::endl;
    std::cout << "F''xx(π, *) = " << F.get_derivative(Derivative::XX)
              << std::endl;
    std::cout << "F''yy(π, *) = " << F.get_derivative(Derivative::YY)
              << std::endl;
    std::cout << "F''xy(π, *) = " << F.get_derivative(Derivative::XY)
              << std::endl;

    std::cout << "G(x, y) = cos(y)" << std::endl;
    std::cout << "G(*, π) = " << G.get_value() << std::endl;
    std::cout << "G'x(*, π) = " << G.get_derivative(Derivative::X) << std::endl;
    std::cout << "G'y(*, π) = " << G.get_derivative(Derivative::Y) << std::endl;
    std::cout << "G''xx(*, π) = " << G.get_derivative(Derivative::XX)
              << std::endl;
    std::cout << "G''yy(*, π) = " << G.get_derivative(Derivative::YY)
              << std::endl;
    std::cout << "G''xy(*, π) = " << G.get_derivative(Derivative::XY)
              << std::endl;

    std::cout << "H(x, y) = exp(y)" << std::endl;
    std::cout << "H(*, π) = " << H.get_value() << std::endl;
    std::cout << "H'x(*, π) = " << H.get_derivative(Derivative::X) << std::endl;
    std::cout << "H'y(*, π) = " << H.get_derivative(Derivative::Y) << std::endl;
    std::cout << "H''xx(*, π) = " << H.get_derivative(Derivative::XX)
              << std::endl;
    std::cout << "H''yy(*, π) = " << H.get_derivative(Derivative::YY)
              << std::endl;
    std::cout << "H''xy(*, π) = " << H.get_derivative(Derivative::XY)
              << std::endl;

    AAD22 A1 = AAD22(Variable::X, M_PI);
    AAD22 A2 = AAD22(Variable::Y, 3);
    //    AAD22 A3 = AAD22(0);
    AAD22 B1 = A1 / (A2 + A1) + A2 * A2;
    AAD22 B2 = A1 - A2 - 1;
    AAD22 B3 = (A1 - A2) * A1 / (A2 * A2);
    //    AAD22 B4 = A1 / A3;
    B1 += 10;

    std::cout << "B1(x, y) = x / (y + x)  + y^2 + 10" << std::endl;
    std::cout << "B1(π, 3) = " << B1.get_value() << std::endl;
    std::cout << "B1'x(π, 3) = " << B1.get_derivative(Derivative::X)
              << std::endl;
    std::cout << "B1'y(π, 3) = " << B1.get_derivative(Derivative::Y)
              << std::endl;
    std::cout << "B1''xx(π, 3) = " << B1.get_derivative(Derivative::XX)
              << std::endl;
    std::cout << "B1''yy(π, 3) = " << B1.get_derivative(Derivative::YY)
              << std::endl;
    std::cout << "B1''xy(π, 3) = " << B1.get_derivative(Derivative::XY)
              << std::endl;

    std::cout << "B2(x, y) = x - y - 1" << std::endl;
    std::cout << "B2(π, 3) = " << B2.get_value() << std::endl;
    std::cout << "B2'x(π, 3) = " << B2.get_derivative(Derivative::X)
              << std::endl;
    std::cout << "B2'y(π, 3) = " << B2.get_derivative(Derivative::Y)
              << std::endl;
    std::cout << "B2''xx(π, 3) = " << B2.get_derivative(Derivative::XX)
              << std::endl;
    std::cout << "B2''yy(π, 3) = " << B2.get_derivative(Derivative::YY)
              << std::endl;
    std::cout << "B2''xy(π, 3) = " << B2.get_derivative(Derivative::XY)
              << std::endl;

    std::cout << "B3(x, y) = (x - y) * x / y^2" << std::endl;
    std::cout << "B3(π, 3) = " << B3.get_value() << std::endl;
    std::cout << "B3'x(π, 3) = " << B3.get_derivative(Derivative::X)
              << std::endl;
    std::cout << "B3'y(π, 3) = " << B3.get_derivative(Derivative::Y)
              << std::endl;
    std::cout << "B3''xx(π, 3) = " << B3.get_derivative(Derivative::XX)
              << std::endl;
    std::cout << "B3''yy(π, 3) = " << B3.get_derivative(Derivative::YY)
              << std::endl;
    std::cout << "B3''xy(π, 3) = " << B3.get_derivative(Derivative::XY)
              << std::endl;

    return 0;
}
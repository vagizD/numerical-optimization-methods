#include "../include/diff.h"

#include <iostream>
#include <cmath>

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
    std::cout << "F''xx(π, *) = " << F.get_derivative(Derivative::XX) << std::endl;
    std::cout << "F''yy(π, *) = " << F.get_derivative(Derivative::YY) << std::endl;
    std::cout << "F''xy(π, *) = " << F.get_derivative(Derivative::XY) << std::endl;

    std::cout << "G(x, y) = cos(y)" << std::endl;
    std::cout << "G(*, π) = " << G.get_value() << std::endl;
    std::cout << "G'x(*, π) = " << G.get_derivative(Derivative::X) << std::endl;
    std::cout << "G'y(*, π) = " << G.get_derivative(Derivative::Y) << std::endl;
    std::cout << "G''xx(*, π) = " << G.get_derivative(Derivative::XX) << std::endl;
    std::cout << "G''yy(*, π) = " << G.get_derivative(Derivative::YY) << std::endl;
    std::cout << "G''xy(*, π) = " << G.get_derivative(Derivative::XY) << std::endl;

    std::cout << "H(x, y) = exp(y)" << std::endl;
    std::cout << "H(*, π) = " << H.get_value() << std::endl;
    std::cout << "H'x(*, π) = " << H.get_derivative(Derivative::X) << std::endl;
    std::cout << "H'y(*, π) = " << H.get_derivative(Derivative::Y) << std::endl;
    std::cout << "H''xx(*, π) = " << H.get_derivative(Derivative::XX) << std::endl;
    std::cout << "H''yy(*, π) = " << H.get_derivative(Derivative::YY) << std::endl;
    std::cout << "H''xy(*, π) = " << H.get_derivative(Derivative::XY) << std::endl;

    return 0;
}
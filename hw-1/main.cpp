#include "exp.hpp"
#include <cmath>
#include <iostream>
#include <cstdint>
#include <iomanip>

int main() {
    float atol_f = 2 * Eps<float>;
    std::cout << std::setprecision(20);

    float range_f[] = {(float)INT32_MIN - 1, -1000, -12.34, -6.5, -1,
                     0,
                     1, 6.5, 12.34, 1000, (float)INT32_MAX + 1};

    for (auto x: range_f) {
        float std_exp = std::exp(x);
        float my_exp =  Exp(x, atol_f);
        float diff = std::abs(std_exp - my_exp);

        if (diff > atol_f) {
            std::cout << "x: " << x << '\n';
            std::cout << "diff: " << diff << '\n';
            std::cout << "exp(x): " << std_exp << '\n';
            std::cout << "my_exp(x): " << my_exp << '\n' << '\n';
        }

//        assert(diff < atol_f);
    }

    double atol_d = 2 * Eps<double>;
    double range_d[] = {(float)INT32_MIN - 1, -1000, -12.34, -6.5, -1,
                     0,
                     1, 6.5, 12.34, 1000, (float)INT32_MAX + 1};

    for (auto x: range_d) {
        double std_exp = std::exp(x);
        double my_exp = Exp(x, atol_d);
        double diff = std::abs(std_exp - my_exp);

        if (diff > atol_d) {
            std::cout << "x: " << x << '\n';
            std::cout << "diff: " << diff << '\n';
            std::cout << "exp(x): " << std_exp << '\n';
            std::cout << "my_exp(x): " << my_exp << '\n' << '\n';
        }

//        assert(diff < atol_d);
    }
}
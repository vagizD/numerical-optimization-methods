#pragma once

#include <array>
#include <cstddef>
#include <iostream>

struct Arg {
    Arg() = default;
    std::array<double, 3> data{};

    explicit Arg(std::array<double, 3> data) : data(data) {
    }

    explicit Arg(double x, double y, double z) : data({x, y, z}) {
    }

    double operator[](size_t idx) const {
        return data[idx];
    }

    double &operator[](size_t idx) {
        return data[idx];
    }

    Arg operator*(double scalar) const {
        double x = data[0] * scalar;
        double y = data[1] * scalar;
        double z = data[2] * scalar;
        return Arg(x, y, z);
    }

    Arg operator*(const Arg &rhs) const {
        double x = data[0] * rhs[0];
        double y = data[1] * rhs[1];
        double z = data[2] * rhs[2];
        return Arg(x, y, z);
    }

    Arg operator+(const Arg &rhs) const {
        double x = data[0] + rhs[0];
        double y = data[1] + rhs[1];
        double z = data[2] + rhs[2];
        return Arg(x, y, z);
    }

    Arg operator-(const Arg &rhs) const {
        double x = data[0] - rhs[0];
        double y = data[1] - rhs[1];
        double z = data[2] - rhs[2];
        return Arg(x, y, z);
    }

    Arg &operator+=(const Arg &rhs) {
        data[0] += rhs[0];
        data[1] += rhs[1];
        data[2] += rhs[2];
        return *this;
    }

    Arg &operator-=(const Arg &rhs) {
        data[0] -= rhs[0];
        data[1] -= rhs[1];
        data[2] -= rhs[2];
        return *this;
    }

    [[nodiscard]] double SquaredNorm() const {
        auto prod = (*this) * (*this);
        return prod[0] + prod[1] + prod[2];
    }

    friend std::ostream &operator<<(std::ostream &os, const Arg arg) {
        for (auto &x : arg.data)
            os << x << ' ';
        os << std::endl;
        return os;
    }
};

class BaseFunction {
public:
    virtual Arg operator()(const double a_t, const Arg &a_y, const Arg &a_dy) {
    }
};

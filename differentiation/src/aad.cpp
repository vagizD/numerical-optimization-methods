#include "aad.h"
#include <cmath>
#include <stdexcept>

// ================ AAD22 GETTERS IMPLEMENTATION ================

double AAD22::get_value() const {
    return m_val;
}

double AAD22::get_derivative(Derivative derivative) const {
    switch (derivative) {
        case Derivative::X:
            return m_d1[0];
        case Derivative::Y:
            return m_d1[1];
        case Derivative::XX:
            return m_d2[0];
        case Derivative::YY:
            return m_d2[1];
        case Derivative::XY:
            return m_d2[2];
        default:
            return 0;
    }
}

// ================ AAD22 OPERATORS IMPLEMENTATION ================

AAD22 AAD22::operator+() const {
    return *this;
}

AAD22 AAD22::operator-() const {
    AAD22 result;
    result.m_val = -m_val;
    result.m_d1 = {-m_d1[0], -m_d1[1]};
    result.m_d2 = {-m_d2[0], -m_d2[1], -m_d2[2]};
    return result;
}

AAD22 &AAD22::operator+=(const AAD22 &rhs) {
    m_val += rhs.m_val;
    m_d1[0] += rhs.m_d1[0];
    m_d1[1] += rhs.m_d1[1];
    m_d2[0] += rhs.m_d2[0];
    m_d2[1] += rhs.m_d2[1];
    m_d2[2] += rhs.m_d2[2];
    return *this;
}

AAD22 AAD22::operator+(const AAD22 &rhs) const {
    AAD22 result = *this;
    result += rhs;
    return result;
}

AAD22 &AAD22::operator-=(const AAD22 &rhs) {
    m_val -= rhs.m_val;
    m_d1[0] -= rhs.m_d1[0];
    m_d1[1] -= rhs.m_d1[1];
    m_d2[0] -= rhs.m_d2[0];
    m_d2[1] -= rhs.m_d2[1];
    m_d2[2] -= rhs.m_d2[2];
    return *this;
}

AAD22 AAD22::operator-(const AAD22 &rhs) const {
    AAD22 result = *this;
    result -= rhs;
    return result;
}

AAD22 &AAD22::operator*=(const AAD22 &rhs) {
    m_d2[0] = m_d2[0] * rhs.m_val + 2 * m_d1[0] * rhs.m_d1[0] + m_val * rhs.m_d2[0];
    m_d2[1] = m_d2[1] * rhs.m_val + 2 * m_d1[1] * rhs.m_d1[1] + m_val * rhs.m_d2[1];
    m_d2[2] = m_d2[2] * rhs.m_val + m_d1[0] * rhs.m_d1[1] + m_d1[1] * rhs.m_d1[0] +
              m_val * rhs.m_d2[2];

    m_d1[0] = m_d1[0] * rhs.m_val + m_val * rhs.m_d1[0];
    m_d1[1] = m_d1[1] * rhs.m_val + m_val * rhs.m_d1[1];

    m_val *= rhs.m_val;

    return *this;
}

AAD22 AAD22::operator*(const AAD22 &rhs) const {
    AAD22 result = *this;
    result *= rhs;
    return result;
}

AAD22 &AAD22::operator/=(const AAD22 &rhs) {
    if (rhs.m_val == 0.0) {
        throw std::runtime_error("Division by zero\n");
    }
    m_d2[0] =
        ((m_d2[0] * rhs.m_val - m_val * rhs.m_d2[0]) * rhs.m_val * rhs.m_val -
         (m_d1[0] * rhs.m_val - m_val * rhs.m_d1[0]) * 2 * rhs.m_val * rhs.m_d1[0]) /
        (rhs.m_val * rhs.m_val * rhs.m_val * rhs.m_val);
    m_d2[1] =
        ((m_d2[1] * rhs.m_val - m_val * rhs.m_d2[1]) * rhs.m_val * rhs.m_val -
         (m_d1[1] * rhs.m_val - m_val * rhs.m_d1[1]) * 2 * rhs.m_val * rhs.m_d1[1]) /
        (rhs.m_val * rhs.m_val * rhs.m_val * rhs.m_val);
    m_d2[2] =
        ((m_d2[2] * rhs.m_val + m_d1[0] * rhs.m_d1[1] - m_d1[1] * rhs.m_d1[0] -
          m_val * rhs.m_d2[2]) *
             rhs.m_val * rhs.m_val -
         (m_d1[0] * rhs.m_val - m_val * rhs.m_d1[0]) * 2 * rhs.m_val * rhs.m_d1[1]) /
        (rhs.m_val * rhs.m_val * rhs.m_val * rhs.m_val);

    m_d1[0] = (m_d1[0] * rhs.m_val - m_val * rhs.m_d1[0]) / (rhs.m_val * rhs.m_val);
    m_d1[1] = (m_d1[1] * rhs.m_val - m_val * rhs.m_d1[1]) / (rhs.m_val * rhs.m_val);

    m_val /= rhs.m_val;

    return *this;
}

AAD22 AAD22::operator/(const AAD22 &rhs) const {
    AAD22 result = *this;
    result /= rhs;
    return result;
}

AAD22 &AAD22::operator+=(const double rhs) {
    m_val += rhs;
    return *this;
}

AAD22 AAD22::operator+(const double rhs) const {
    AAD22 result = *this;
    result += rhs;
    return result;
}

AAD22 &AAD22::operator-=(const double rhs) {
    m_val -= rhs;
    return *this;
}

AAD22 AAD22::operator-(const double rhs) const {
    AAD22 result = *this;
    result -= rhs;
    return result;
}

AAD22 &AAD22::operator*=(const double rhs) {
    m_val *= rhs;
    m_d1[0] *= rhs;
    m_d1[1] *= rhs;
    m_d2[0] *= rhs;
    m_d2[1] *= rhs;
    m_d2[2] *= rhs;
    return *this;
}

AAD22 AAD22::operator*(const double rhs) const {
    AAD22 result = *this;
    result *= rhs;
    return result;
}

AAD22 &AAD22::operator/=(const double rhs) {
    if (rhs == 0.0) {
        throw std::runtime_error("Division by zero\n");
    }
    m_val /= rhs;
    m_d1[0] /= rhs;
    m_d1[1] /= rhs;
    m_d2[0] /= rhs;
    m_d2[1] /= rhs;
    m_d2[2] /= rhs;
    return *this;
}

AAD22 AAD22::operator/(const double rhs) const {
    AAD22 result = *this;
    result /= rhs;
    return result;
}

// ================ AAD22 FUNCTIONS IMPLEMENTATION ================

AAD22 sin(const AAD22 &arg) {
    AAD22 res;
    double arg_cos = std::cos(arg.m_val);
    double arg_sin = std::sin(arg.m_val);
    res.m_val = arg_sin;
    res.m_d1[0] = arg_cos * arg.m_d1[0];
    res.m_d1[1] = arg_cos * arg.m_d1[1];
    res.m_d2[0] = arg_cos * arg.m_d2[0] - arg_sin * arg.m_d1[0] * arg.m_d1[0];
    res.m_d2[1] = arg_cos * arg.m_d2[1] - arg_sin * arg.m_d1[1] * arg.m_d1[1];
    res.m_d2[2] = arg_cos * arg.m_d2[2] - arg_sin * arg.m_d1[0] * arg.m_d1[1];
    return res;
}

AAD22 cos(const AAD22 &arg) {
    AAD22 res;
    double arg_cos = std::cos(arg.m_val);
    double arg_sin = std::sin(arg.m_val);
    res.m_val = arg_cos;
    res.m_d1[0] = -arg_sin * arg.m_d1[0];
    res.m_d1[1] = -arg_sin * arg.m_d1[1];
    res.m_d2[0] = -arg_sin * arg.m_d2[0] - arg_cos * arg.m_d1[0] * arg.m_d1[0];
    res.m_d2[1] = -arg_sin * arg.m_d2[1] - arg_cos * arg.m_d1[1] * arg.m_d1[1];
    res.m_d2[2] = -arg_sin * arg.m_d2[2] - arg_cos * arg.m_d1[0] * arg.m_d1[1];
    return res;
}

AAD22 exp(const AAD22 &arg) {
    AAD22 res;
    double arg_exp = std::exp(arg.m_val);
    res.m_val = arg_exp;
    res.m_d1[0] = arg_exp * arg.m_d1[0];
    res.m_d1[1] = arg_exp * arg.m_d1[1];
    res.m_d2[0] = arg_exp * (arg.m_d2[0] + arg.m_d1[0] * arg.m_d1[0]);
    res.m_d2[1] = arg_exp * (arg.m_d2[1] + arg.m_d1[1] * arg.m_d1[1]);
    res.m_d2[2] = arg_exp * (arg.m_d2[2] + arg.m_d1[0] * arg.m_d1[1]);
    return res;
}

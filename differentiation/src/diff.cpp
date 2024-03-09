#include "../include/diff.h"

#include <cmath>




// ================ AAD22 GETTERS IMPLEMENTATION ================

double AAD22::get_value() const {
    return m_val;
}

double AAD22::get_derivative(Derivative derivative) const {
    switch(derivative) {
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

// ================ AAD22 FUNCTIONS IMPLEMENTATION ================

AAD22 sin(const AAD22& arg) {
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

AAD22 cos(const AAD22& arg) {
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

AAD22 exp(const AAD22& arg) {
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
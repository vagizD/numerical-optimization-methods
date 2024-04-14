#pragma once

#include <array>

constexpr double g = 9.806555;
constexpr double R = 287.0528;

constexpr std::array<double, 4> r = {6.5e-3, 0, -1e-3, -2.8e-3};
constexpr std::array<double, 5> h = {0, 11e3, 20e3, 32e3, 47e3};

// We use coefficients for rkf45, formula 2 table III in Fehlberg
constexpr std::array<double, 1> k2_consts = { 1./4. };
constexpr std::array<double, 2> k3_consts = { 3./32., 9./32. };
constexpr std::array<double, 3> k4_consts = { 1932./2197., -7200./2197., 7296./2197. };
constexpr std::array<double, 4> k5_consts = { 439./216., -8., 3680./513., -845./4104. };
constexpr std::array<double, 5> k6_consts = { -8./27., 2., -3544./2565., 1859./4104., -11./40. };
constexpr std::array<double, 6> delta_consts = { -1./360., 0., 128./4275., 2197./75240., -1./50., -2./55.};
constexpr std::array<double, 6> gamma_consts = { 16./135, 0., 6656./12825., 28561./56430., -9./50., 2./55.};
#include "sam.h"
#include "constants.h"

#include <cassert>
#include <cmath>

Pressure::Pressure(double p0, double t0) {
    t[0] = t0, p[0] = p0;
    for (int i = 1; i < 4; ++i) {
        t[i] = t[i - 1] - r[i - 1] * (h[i] - h[i - 1]);
        p[i] = p[i - 1] * std::exp(g / (R * r[i - 1]) * std::log(1 - r[i - 1] * (h[i] - h[i - 1]) / t[i - 1]));
    }
}

int getLayer(double height) {
    assert(0 <= height && height <= h.back());
    if (height <= h[1]) {
        return 0;
    } else if (height <= h[2]) {
        return 1;
    } else if (height <= h[3]) {
        return 2;
    }
    return 3;
}

double Pressure::operator()(double height) {
    int l = getLayer(height);
    if (r[l] != 0) {
        return p[l] * std::exp(g / (R * r[l]) * std::log(1 - r[l] * (height - h[l]) / t[l]));
    }
    return p[l] * std::exp(g * (h[l] - height) / (R * t[l]));
}

double Density::operator()(double height) {
    int l = getLayer(height);
    return pressure(height) / (R * (pressure.t[l] - r[l] * (height - h[l])));
}

double AirSpeed::operator()(double height) {
    return std::sqrt(density.pressure(height) / density(height));
}











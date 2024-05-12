#pragma once

#include <array>
#include <fstream>
#include <string>
#include <vector>
#include "function.h"

class Point3d {
public:
    Point3d(double a_x, double a_y, double a_z, double a_t)
        : m_x(a_x), m_y(a_y), m_z(a_z), m_t(a_t){};

    [[nodiscard]] double get_x() const {
        return m_x;
    }

    [[nodiscard]] double get_y() const {
        return m_y;
    }

    [[nodiscard]] double get_z() const {
        return m_z;
    }

    [[nodiscard]] double get_t() const {
        return m_t;
    }

private:
    double m_x;
    double m_y;
    double m_z;
    double m_t;
};

class SatelliteObserver {
private:
    std::vector<Point3d> m_projectile;
    double m_x0;
    double m_y0;
    double m_z0;
    double m_t0 = -1;
    double m_tEnd = -1;
    double m_tRes;
    double m_orbit;
    double m_R;
    bool m_save;

public:
    SatelliteObserver(
        double a_x0,
        double a_y0,
        double a_z0,
        double a_orbit,
        double a_R,
        double a_tRes = -1
    )
        : m_x0(a_x0),
          m_y0(a_y0),
          m_z0(a_z0),
          m_orbit(a_orbit),
          m_R(a_R),
          m_tRes(a_tRes){};
    bool make_decision(double a_cur_t, Arg a_cur_y, Arg a_cur_dt, bool verbose);
    void save_projectile(std::string &filename);
};
#pragma once

#include <array>
#include <fstream>
#include <string>
#include <vector>

class Point {
public:
    Point(double a_x, double a_y, double a_t) : m_x(a_x), m_y(a_y), m_t(a_t){};

    double get_x() const {
        return m_x;
    }

    double get_y() const {
        return m_y;
    }

    double get_t() const {
        return m_t;
    }

private:
    double m_x;
    double m_y;
    double m_t;
};

class SimpleObserver {
private:
    std::vector<Point> m_projectile;
    double m_x0;
    double m_y0;
    double m_t0 = -1;
    double m_tEnd = -1;
    bool m_save;

public:
    double m_y_max = -1;
    SimpleObserver(double a_x0, double a_y0) : m_x0(a_x0), m_y0(a_y0){};
    bool make_decision(  // observation + stop criteria
        double a_cur_t,
        std::array<double, 4> a_cur_y,
        bool verbose
    );
    void save_projectile(std::string &filename);
};

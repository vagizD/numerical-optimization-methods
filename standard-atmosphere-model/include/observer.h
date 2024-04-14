#pragma once

#include <fstream>
#include <string>
#include <vector>

class Point {
public:
    Point(int a_x, int a_y, int a_t) : m_x(a_x), m_y(a_y), m_t(a_t){};

    double get_x() {
        return m_x;
    }

    double get_y() {
        return m_y;
    }

    double get_t() {
        return m_t;
    }

private:
    int m_x;
    int m_y;
    int m_t;
};

class SimpleObserver {
private:
    std::vector<Point> m_projectile;
    double m_x0;
    double m_y0;
    double m_t0 = -1;
    double m_tEnd = -1;

public:
    SimpleObserver(double a_x0, double a_y0) : m_x0(a_x0), m_y0(a_y0){};
    // a_n -- RHS::N, a_cur_y -- vector y
    bool operator()(
        double a_cur_t,
        int a_n,
        const double *a_cur_y,
        bool verbose
    );  // observation + stop criteria

    void save_projectile(const std::string &filename) {
        std::ofstream fhandle;
        fhandle.open(filename);
        fhandle << "t, x, y\n";
        for (auto &p : m_projectile) {
            fhandle << p.get_t() << ", " << p.get_x() << ", " << p.get_y() << '\n';
        }
        fhandle.close();
    }
};

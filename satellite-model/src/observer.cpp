#include "observer.h"
#include <iostream>

bool SatelliteObserver::make_decision(  // observation + stop criteria
    double a_cur_t,
    Arg a_cur_y,
    Arg a_cur_dy,
    bool verbose
) {
    if (a_cur_y.SquaredNorm() < m_R) {
        std::cout << "CRUSHED!";
        return false;
    }
    if (std::abs(a_cur_y[0]) > m_orbit or std::abs(a_cur_y[1]) > m_orbit or
        std::abs(a_cur_y[2]) > m_orbit) {
        std::cout << "LEFT ORBIT";
        return false;
    }
    if (a_cur_t > m_tRes) {
        return false;
    }

    if (m_t0 == -1)
        m_t0 = a_cur_t;

    m_tEnd = a_cur_t;
    m_projectile.emplace_back(a_cur_y[0], a_cur_y[1], a_cur_y[2], a_cur_t);

    if (verbose) {
        std::cout << "State: t = " << a_cur_t << ", x = " << a_cur_y[0]
                  << ", y = " << a_cur_y[1] << ", z = " << a_cur_y[2]
                  << ", R = " << a_cur_y.SquaredNorm() << ", V_x = " << a_cur_dy[0]
                  << ", V_y = " << a_cur_dy[1] << ", V_z = " << a_cur_dy[2] << std::endl;
    }
    return true;
}

void SatelliteObserver::save_projectile(std::string &filename) {
    std::ofstream fhandle;
    fhandle.open(filename);
    fhandle << "t,x,y,z\n";
    for (auto &p : m_projectile) {
        fhandle << p.get_t() << ", " << p.get_x() << ", " << p.get_y() << ", "
                << p.get_z() << '\n';
    }
    fhandle.close();
}
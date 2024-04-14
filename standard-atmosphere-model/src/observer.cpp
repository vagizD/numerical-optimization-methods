#include "observer.h"
#include <array>
#include <iostream>

bool SimpleObserver::make_decision(  // observation + stop criteria
    double a_cur_t,
    std::array<double, 4> a_cur_y,
    bool verbose
) {
    if (a_cur_y[2] < 0) {  // stop if y coordinate <= 0
        return false;
    }

    if (m_t0 == -1)
        m_t0 = a_cur_t;
    m_tEnd = a_cur_t;
    m_projectile.emplace_back(a_cur_y[0], a_cur_y[2], a_cur_t);

    if (verbose) {
        std::cout << "State: t = " << a_cur_t << ", x = " << a_cur_y[0]
                  << ", y = " << a_cur_y[2] << std::endl;
    }
    return true;
}

void SimpleObserver::save_projectile(std::string &filename) {
    std::ofstream fhandle;
    fhandle.open(filename);
    fhandle << "t,x,y\n";
    for (auto &p : m_projectile) {
        fhandle << p.get_t() << ", " << p.get_x() << ", " << p.get_y() << '\n';
    }
    fhandle.close();
}

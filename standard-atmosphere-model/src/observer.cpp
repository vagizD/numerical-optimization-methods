#include "observer.h"


bool SimpleObserver::operator()(double a_cur_t, int a_n, const double *const a_cur_y) {
    if (m_t0 == -1) m_t0 = a_cur_t;
    m_tEnd = a_cur_t;
    m_projectile.emplace_back(a_cur_y[0], a_cur_y[2], a_cur_t);

    return a_cur_y[2] > 0;  // stop if y coordinate <= 0
}

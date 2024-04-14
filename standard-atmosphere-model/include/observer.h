#pragma once

class SimpleObserver {
public:
    bool operator()(
        double a_cur_t,
        int a_n,
        const double *const a_cur_y
    ) {                               // a_n -- RHS::N, a_cur_y -- vector y
        return a_cur_y[a_n - 2] > 0;  // stop if y coordinate <= 0
    }
};
#pragma once

#include <set>
#include "algebra.hpp"
#include "constants.hpp"
#include "sturm.hpp"
#include "system.hpp"


namespace ADAAI {

template<typename F, size_t N, size_t M>
std::set<Zero> backtrack(const System<F, N, M>& S, const size_t curIdx, Zero zero) {
    if (curIdx == 0) {
        return {zero};
    }

    UnivariatePoly<F, M> g = S[curIdx].substitute(curIdx, zero);

    SturmF V = buildSturmSequence(g);
    std::vector<HalfInterval<F>> singleZeroHalfIntervals = findHalfIntervals(V, INITIAL_STEP);

    std::set<Zero> newZeros;
    for (auto& halfInterval: singleZeroHalfIntervals) {
        shrinkHalfInterval(halfInterval, INTERVAL_SIZE);
        const double x = approximateZero(halfInterval);
        Zero curZero = zero + Zero{x};
        std::set<Zero> foundZeros = backtrack(S, curIdx - 1, curZero);
        newZeros += foundZeros;
    }
    return newZeros;
}

}  // namespace ADAAI


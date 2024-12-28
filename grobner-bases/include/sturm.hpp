#pragma once

#include "polynomial.hpp"

namespace ADAAI {

template<typename F>
class HalfInterval {
public:
    F a;
    F b;

    explicit HalfInterval(F a, F b): a(a), b(b) {}
};

template<typename F, size_t N, size_t M>
class SturmF {
public:
    std::vector<MultivariatePoly<F, N, M>> polynomials;

    explicit SturmF(const std::vector<MultivariatePoly<F, N, M>>& polynomials):
    polynomials(polynomials) {}
};

}  // namespace ADAAI

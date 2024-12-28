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

template<typename F, size_t N>
class SturmF {
public:
    std::vector<UnivariatePoly<F, N>> polynomials;

    explicit SturmF(const std::vector<UnivariatePoly<F, N>>& polynomials):
    polynomials(polynomials) {}
};

}  // namespace ADAAI

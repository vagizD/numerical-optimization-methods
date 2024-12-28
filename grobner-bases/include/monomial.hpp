#pragma once

#include <array>

namespace ADAAI {

template <size_t N>
class Monomial {
public:
    std::array<size_t, N> degrees;  // degree of x_1, ..., x_N

    explicit Monomial(std::array<size_t, N>& degrees):
    degrees(degrees) {}
};

template <typename F, size_t N>
class Term {
public:
    F coefficient;
    Monomial<N> monomial;

    explicit Term(F coefficient, Monomial<N>& monomial):
    coefficient(coefficient), monomial(monomial) {}

    explicit Term(F coefficient, std::array<size_t, N>& degrees):
    coefficient(coefficient), monomial(degrees) {}
};

template<typename F, size_t N>
typedef bool termOrdering(const Term<F, N>& t1, const Term<F, N>& t2);

}  // namespace ADAAI

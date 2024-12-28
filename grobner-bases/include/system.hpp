#pragma once

#include <vector>
#include "polynomial.hpp"

namespace ADAAI {

template<typename F, size_t N, size_t M>
class System {
public:
    std::vector<MultivariatePoly<F, N, M>> polynomials;

    explicit System(const std::vector<MultivariatePoly<F, N, M>>& polynomials):
    polynomials(polynomials) {}

    void& operator[](const size_t idx) {
        return polynomials[idx];
    }
};

}  // namespace ADAAI

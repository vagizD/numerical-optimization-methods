#pragma once

#include <algorithm>
#include <array>

namespace ADAAI {
template <typename F, size_t N>
class Poly {
private:
    std::array<F, N> coefficients;

public:
    explicit Poly(const std::array<F, N> &a_coefficients)
        : coefficients(a_coefficients) {
    }

    F eval(F a_x) const {
        F v = 0.0;
        F p = 1.0;
        std::vector<F> terms;
        // evaluate poly(a_x)
        for (F coef : coefficients) {
            terms.push_back(p * coef);
            p *= a_x;
        }
        // summate ascending modulus of the terms
        std::reverse(terms.begin(), terms.end());
        for (F term : terms)
            v += term;
        return v;
    }
};
}  // namespace ADAAI

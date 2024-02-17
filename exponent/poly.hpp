#pragma once

#include <algorithm>
#include <array>
#include <cstdio>
#include <vector>

namespace ADAAI {
template <typename F, size_t N>
class Poly {
private:
    std::array<F, N> coefficients;
    size_t deg = 0;

public:
    explicit Poly(const std::array<F, N> &a_coefficients, size_t deg)
        : coefficients(a_coefficients), deg(deg) {
    }

    [[nodiscard]] F eval(F a_x) const {
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

    Poly<F, N> &operator+=(Poly<F, N> &other) {
        for (size_t i = 0; i < std::max(other.deg, deg))
            coefficients[i] += other.coefficients[i];
        size = std::max(other.deg, deg);
        return *this;
    };

    Poly<F, N> &operator*=(Poly<F, N> &other) {
        std::array<F, N> mul;
        static_assert(deg + other.deg < N);
        for (size_t i = 0; i <= deg; i++) {
            for (size_t j = 0; j <= other.deg; j++) {
                mul[i + j] = coefficients[i] * other.coefficients[j];
            }
        }
        deg += other.deg;
        coefficients = std::move(mul);
        return *this;
    };

    // => Horner's method <=
    std::pair<Poly<F, N>, Poly<F, N>> Poly Horner(const Poly<F, N> &other) {
        std::array<F, N> res = coefficients;
        int cur = deg;
        while (cur - other.deg >= 0) {
            res[cur] = coefficients[cur] / other.coefficients[other.deg];
            for (int i = 1; i <= other.deg; ++i) {
                res[cur - i] -= res[cur] * other.coefficients[other.deg - i];
            }
            cur--;
        }

        int remDeg = other.deg - 1;
        while (res[remDeg] == 0 && remDeg >= 0)
            remDeg--;

        std::array<F, N> div;
        for (int i = other.deg; i <= deg; ++i) {
            div[i - other.deg] = res[i];
        }

        std::array<F, N> rem;
        for (int i = 0; i <= remDeg; ++i) {
            rem[i] = res[i];
        }

        return std::make_pair(Poly(div, deg - other.deg), Poly(rem, remDeg));
    }

    friend std::ostream &operator<<(std::ostream &os, const Poly<F, N> p) {
        for (auto &x : p.coefficients) {
            os << x << ' ';
        }

        return os;
    }
};
}  // namespace ADAAI
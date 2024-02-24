#pragma once

#include <algorithm>
#include <array>
#include <cassert>
#include <cstdio>
#include <vector>
#include "constants.hpp"

namespace ADAAI {
template <typename F, size_t N>
class Poly {
public:
    std::array<F, N> coefficients;
    size_t deg;

    explicit Poly(const std::array<F, N> &a_coefficients)
        : coefficients(a_coefficients) {
        deg = N - 1;
        while (std::abs(a_coefficients[deg]) < (F)10.0 * ADAAI::Eps<F> &&
               deg != 0)
            coefficients[deg--] = 0;
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

    Poly<F, N> operator-(const Poly<F, N> &other) {
        std::array<F, N> sub = coefficients;
        for (size_t i = 0; i <= std::max(other.deg, deg); i++)
            sub[i] -= other.coefficients[i];
        return Poly<F, N>(sub);
    };

    Poly<F, N> operator+(const Poly<F, N> &other) {
        std::array<F, N> sum = coefficients;
        for (size_t i = 0; i <= std::max(other.deg, deg); i++)
            sum[i] += other.coefficients[i];
        return Poly<F, N>(sum);
    };

    Poly<F, N> &operator+=(const Poly<F, N> &other) {
        for (size_t i = 0; i <= std::max(other.deg, deg); i++)
            coefficients[i] += other.coefficients[i];
        deg = std::max(other.deg, deg);
        return *this;
    };

    Poly<F, N> operator*(const Poly<F, N> &other) {
        std::array<F, N> mul = {};
        assert(deg + other.deg < N);
        for (size_t i = 0; i <= deg; i++) {
            for (size_t j = 0; j <= other.deg; j++) {
                mul[i + j] += coefficients[i] * other.coefficients[j];
            }
        }
        return Poly<F, N>(mul);
    };

    Poly<F, N> &operator*=(const Poly<F, N> &other) {
        std::array<F, N> mul = {};
        assert(deg + other.deg < N);
        for (size_t i = 0; i <= deg; i++) {
            for (size_t j = 0; j <= other.deg; j++) {
                mul[i + j] += coefficients[i] * other.coefficients[j];
            }
        }
        deg += other.deg;
        coefficients = std::move(mul);
        return *this;
    };

    // => Horner's method <=
    Poly<F, N> operator/(const Poly<F, N> &other) {
        std::array<F, N> res = coefficients;
        int cur = deg;
        while (cur >= static_cast<int>(other.deg)) {
            res[cur] = res[cur] / other.coefficients[other.deg];
            for (int i = 1; i <= other.deg; ++i) {
                res[cur - i] -= res[cur] * other.coefficients[other.deg - i];
            }
            cur--;
        }

        std::array<F, N> div = {};
        for (int i = other.deg; i <= deg; ++i) {
            div[i - other.deg] = res[i];
        }

        return Poly<F, N>(div);
    }

    friend std::ostream &operator<<(std::ostream &os, const Poly<F, N> p) {
        for (auto &x : p.coefficients) {
            os << x << ' ';
        }

        return os;
    }
};
}  // namespace ADAAI
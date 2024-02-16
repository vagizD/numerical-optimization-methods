#pragma once

#include <algorithm>
#include <array>
#include <cstdio>

namespace ADAAI {
template <typename F, size_t N>
class Poly {
private:
    std::array<F, N> coefficients;
    size_t size = 0;

public:
    explicit Poly(const std::array<F, N> &a_coefficients, size_t s)
        : coefficients(a_coefficients), size(s) {
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

    Poly<F, N> &operator+=(Poly<F, N> &poly1) {
        for (size_t i = 0; i < N; i++) {
            if (i >= std::max(poly1.size, this->size)) {
                this->coefficients[i] = 0;
            } else {
                this->coefficients[i] += poly1.coefficients[i];
            }
        }

        this->size = std::max(poly1.size, this->size);

        return *this;
    };

    Poly<F, N> &operator*=(Poly<F, N> &poly1) {
        Poly<F, N> copy_t(this->coefficients, this->size);

        for (int i = 0; i < this->size; i++) {
            this->coefficients[i] = 0;
        }

        if (this->size == 0 || poly1.size == 0) {
            this->size = 0;
            return *this;
        }

        for (size_t i = 0; i < this->size; i++) {
            for (size_t j = 0; j < poly1.size; j++) {
                this->coefficients[i + j] +=
                    copy_t.coefficients[i] * poly1.coefficients[j];
                //                    if(this->coefficients[i+j] <= 10e-7){
                //                        this->coefficients[i+j] = 0;
                //                    }
            }
        }

        this->size = this->size + poly1.size - 1;

        return *this;
    };

    friend std::ostream &operator<<(std::ostream &os, const Poly<F, N> p) {
        for (auto &x : p.coefficients) {
            os << x << ' ';
        }

        return os;
    }
};
}  // namespace ADAAI
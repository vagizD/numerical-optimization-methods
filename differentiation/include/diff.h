#pragma once

#include <array>

enum class DiffMethod { Stencil3, Stencil3Extra, Stencil5, Stencil5Extra, FwdADD };

enum class Derivative { X, Y, XX, YY, XY };

enum class Variable { X, Y };

template <Derivative D, DiffMethod M, typename Callable>
double Differentiator(Callable F, double x, double y);

class AAD22 {
   public:
    AAD22() : m_val(0){};
    explicit AAD22(double v) : m_val(v) {}
    AAD22(Variable var, double v) : m_val(v) {
        if (var == Variable::X) {
            m_d1[0] = 1;
        } else if (var == Variable::Y) {
            m_d1[1] = 1;
        }
    }

    AAD22 operator+() const;
    AAD22 operator-() const;

    AAD22& operator+=(const AAD22& rhs);
    AAD22& operator-=(const AAD22& rhs);
    AAD22& operator*=(const AAD22& rhs);
    AAD22& operator/=(const AAD22& rhs);

    AAD22 operator+(const AAD22& rhs) const;
    AAD22 operator-(const AAD22& rhs) const;
    AAD22 operator*(const AAD22& rhs) const;
    AAD22 operator/(const AAD22& rhs) const;

    AAD22& operator+=(double rhs);
    AAD22& operator-=(double rhs);
    AAD22& operator*=(double rhs);
    AAD22& operator/=(double rhs);

    AAD22 operator+(double rhs) const;
    AAD22 operator-(double rhs) const;
    AAD22 operator*(double rhs) const;
    AAD22 operator/(double rhs) const;

    friend AAD22 sin(const AAD22& arg);
    friend AAD22 cos(const AAD22& arg);
    friend AAD22 exp(const AAD22& arg);

    [[nodiscard]] double get_value() const;
    [[nodiscard]] double get_derivative(Derivative derivative) const;

   private:
    double m_val;
    std::array<double, 2> m_d1 = {};
    std::array<double, 3> m_d2 = {};
};

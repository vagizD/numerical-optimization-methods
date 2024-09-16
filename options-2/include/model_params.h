#pragma once
#include <cassert>

class ModelParams {
public:
    double r;        // risk-interest rate
    double d;        // dividend rate
    double sigma;    // volatility coefficient
    double kappa;    // stochastic volatility coefficient
    double beta;     // constant elasticity of variance
    double epsilon;  // stochastic volatility coefficient
    double rho;      // correlation of brownian differentials
    double theta;    //

    constexpr ModelParams(
        const double r,
        const double d,
        const double sigma,
        const double kappa,
        const double beta,
        const double epsilon,
        const double rho,
        const double theta
    )
        : r(r),
          d(d),
          sigma(sigma),
          kappa(kappa),
          beta(beta),
          epsilon(epsilon),
          rho(rho),
          theta(theta) {
        assert((0 <= rho) && (rho < 1));
        assert(kappa > 0);
        assert(theta > 0);
        assert(epsilon > 0);
        assert(sigma > 0);
        assert((0 < beta) && (beta <= 1));
    };
};

#include <iomanip>
#include <iostream>
#include "constants.hpp"
#include "exp.hpp"
#include "test.hpp"

int main() {
    std::cout << std::setprecision(10);

    const size_t N = 50;  // Capacity

    // we test over [-Ln2/2, Ln2/2] in steps of 0.001
    {
        float step = 0.01;
        float r = ADAAI::Ln2<float>() / 2;
        std::cout << "===RESULTS FOR FLOAT===" << std::endl;
        std::cout << "=> EPS is set to be "
                  << 10.0 * ADAAI::Eps<float> << std::endl;
        {
            auto [absError, relError] = ADAAI::makeTests<
                float>(-r, r, step, ADAAI::Exp<ADAAI::MethodE::Pade, float, N>);
            std::cout << "=> PADE   | Max absolute error: " << absError
                      << std::endl;
            std::cout << "=> PADE   | Max relative error: " << relError
                      << std::endl;
        }
        {
            auto [absError, relError] = ADAAI::makeTests<
                float>(-r, r, step, ADAAI::Exp<ADAAI::MethodE::Taylor, float, N>);
            std::cout << "=> Taylor | Max absolute error: " << absError
                      << std::endl;
            std::cout << "=> Taylor | Max relative error: " << relError
                      << std::endl;
        }
        std::cout << std::endl;
    }
    {
        double step = 0.01;
        double r = ADAAI::Ln2<double>() / 2;
        std::cout << "===RESULTS FOR DOUBLE===" << std::endl;
        std::cout << "=> EPS is set to be "
                  << 10.0 * ADAAI::Eps<double> << std::endl;
        ;
        {
            auto [absError, relError] = ADAAI::makeTests<
                double>(-r, r, step, ADAAI::Exp<ADAAI::MethodE::Pade, double, N>);
            std::cout << "=> PADE   | Max absolute error: " << absError
                      << std::endl;
            std::cout << "=> PADE   | Max relative error: " << relError
                      << std::endl;
        }
        {
            auto [absError, relError] = ADAAI::makeTests<
                double>(-r, r, step, ADAAI::Exp<ADAAI::MethodE::Taylor, double, N>);
            std::cout << "=> Taylor | Max absolute error: " << absError
                      << std::endl;
            std::cout << "=> Taylor | Max relative error: " << relError
                      << std::endl;
        }
        std::cout << std::endl;
    }
    {
        long double step = 0.01;
        long double r = ADAAI::Ln2<long double>() / 2;
        std::cout << "===RESULTS FOR LONG DOUBLE===" << std::endl;
        std::cout << "=> EPS is set to be "
                  << 10.0 * ADAAI::Eps<long double> << std::endl;
        {
            auto [absError, relError] = ADAAI::makeTests<
                long double>(-r, r, step, ADAAI::Exp<ADAAI::MethodE::Pade, long double, N>);
            std::cout << "=> PADE   | Max absolute error: " << absError
                      << std::endl;
            std::cout << "=> PADE   | Max relative error: " << relError
                      << std::endl;
        }
        {
            auto [absError, relError] = ADAAI::makeTests<
                long double>(-r, r, step, ADAAI::Exp<ADAAI::MethodE::Taylor, long double, N>);
            std::cout << "=> Taylor | Max absolute error: " << absError
                      << std::endl;
            std::cout << "=> Taylor | Max relative error: " << relError
                      << std::endl;
        }
    }

    // LOCAL RESULTS
    // ==============================================
    // ===RESULTS FOR FLOAT===
    // => EPS is set to be 1.192092896e-06
    // => PADE   | Max absolute error: 1.192092896e-07
    // => PADE   | Max relative error: 1.25017678e-07
    // => Taylor | Max absolute error: 8.344650269e-07
    // => Taylor | Max relative error: 1.115842679e-06
    //
    // ===RESULTS FOR DOUBLE===
    // => EPS is set to be 2.220446049e-15
    // => PADE   | Max absolute error: 4.440892099e-16
    // => PADE   | Max relative error: 3.259027603e-16
    // => Taylor | Max absolute error: 1.554312234e-15
    // => Taylor | Max relative error: 2.11616766e-15
    //
    // ===RESULTS FOR LONG DOUBLE===
    // => EPS is set to be 1.084202172e-18
    // => PADE   | Max absolute error: 8.67361738e-19
    // => PADE   | Max relative error: 6.134076417e-19
    // => Taylor | Max absolute error: 7.589415207e-19
    // => Taylor | Max relative error: 9.575618792e-19

    std::cout << "\n=> Tests on ranges larger than [-10, 10] <=\n\n";
    // type: float
    // range: [-20, 20]
    // step: 0.001
    {
        float step = 0.01;
        float r = 20;
        std::cout << "===RESULTS FOR FLOAT===" << std::endl;
        std::cout << "=> EPS is set to be "
                  << 10.0 * ADAAI::Eps<float> << std::endl;
        {
            auto [absError, relError] = ADAAI::makeTests<
                float>(-r, r, step, ADAAI::Exp<ADAAI::MethodE::Pade, float, N>);
            std::cout << "=> PADE   | Max absolute error: " << absError
                      << std::endl;
            std::cout << "=> PADE   | Max relative error: " << relError
                      << std::endl;
        }
        {
            auto [absError, relError] = ADAAI::makeTests<
                float>(-r, r, step, ADAAI::Exp<ADAAI::MethodE::Taylor, float, N>);
            std::cout << "=> Taylor | Max absolute error: " << absError
                      << std::endl;
            std::cout << "=> Taylor | Max relative error: " << relError
                      << std::endl;
        }
        std::cout << std::endl;
    }
    // type: double
    // range: [-40, 40]
    // step: 0.001
    {
        double step = 0.01;
        double r = 40;
        std::cout << "===RESULTS FOR DOUBLE===" << std::endl;
        std::cout << "=> EPS is set to be "
                  << 10.0 * ADAAI::Eps<double> << std::endl;
        ;
        {
            auto [absError, relError] = ADAAI::makeTests<
                double>(-r, r, step, ADAAI::Exp<ADAAI::MethodE::Pade, double, N>);
            std::cout << "=> PADE   | Max absolute error: " << absError
                      << std::endl;
            std::cout << "=> PADE   | Max relative error: " << relError
                      << std::endl;
        }
        {
            auto [absError, relError] = ADAAI::makeTests<
                double>(-r, r, step, ADAAI::Exp<ADAAI::MethodE::Taylor, double, N>);
            std::cout << "=> Taylor | Max absolute error: " << absError
                      << std::endl;
            std::cout << "=> Taylor | Max relative error: " << relError
                      << std::endl;
        }
        std::cout << std::endl;
    }
    // type: long double
    // range: [-80, 80]
    // step: 0.001
    {
        long double step = 0.01;
        long double r = 80;
        std::cout << "===RESULTS FOR LONG DOUBLE===" << std::endl;
        std::cout << "=> EPS is set to be "
                  << 10.0 * ADAAI::Eps<long double> << std::endl;
        {
            auto [absError, relError] = ADAAI::makeTests<
                long double>(-r, r, step, ADAAI::Exp<ADAAI::MethodE::Pade, long double, N>);
            std::cout << "=> PADE   | Max absolute error: " << absError
                      << std::endl;
            std::cout << "=> PADE   | Max relative error: " << relError
                      << std::endl;
        }
        {
            auto [absError, relError] = ADAAI::makeTests<
                long double>(-r, r, step, ADAAI::Exp<ADAAI::MethodE::Taylor, long double, N>);
            std::cout << "=> Taylor | Max absolute error: " << absError
                      << std::endl;
            std::cout << "=> Taylor | Max relative error: " << relError
                      << std::endl;
        }
    }

    return 0;
}

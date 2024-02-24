#include <iomanip>
#include <iostream>
#include "constants.hpp"
#include "test.hpp"

int main() {
    std::cout << std::setprecision(10);

    // we test over [-Ln2/2, Ln2/2] in steps of 0.001
    {
        float step = 0.01;
        float r = ADAAI::Ln2<float>() / 2;
        ADAAI::runTests<float>(-r, r, step, "FLOAT");
    }
    {
        double step = 0.01;
        double r = ADAAI::Ln2<double>() / 2;
        ADAAI::runTests<double>(-r, r, step, "DOUBLE");
    }
    {
        long double step = 0.01;
        long double r = ADAAI::Ln2<long double>() / 2;
        ADAAI::runTests<long double>(-r, r, step, "LONG DOUBLE");
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
        ADAAI::runTests<float>(-r, r, step, "FLOAT");
    }
    // type: double
    // range: [-40, 40]
    // step: 0.001
    {
        double step = 0.01;
        double r = 40;
        ADAAI::runTests<double>(-r, r, step, "DOUBLE");
    }
    // type: long double
    // range: [-80, 80]
    // step: 0.001
    {
        long double step = 0.01;
        long double r = 80;
        ADAAI::runTests<long double>(-r, r, step, "LONG DOUBLE");
    }

    return 0;
}

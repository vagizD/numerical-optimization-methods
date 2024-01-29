#include "exp.hpp"
#include <cmath>
#include <iostream>
#include <iomanip>

int main() {

    // Testing for float values

    // Testing for x <= 0
    float negFloatMaxErr = 0.0f; // absolute error
    float negFloatTests[] = {
            -20.        , -19.59183673, -19.18367347, -18.7755102 ,
            -18.36734694, -17.95918367, -17.55102041, -17.14285714,
            -16.73469388, -16.32653061, -15.91836735, -15.51020408,
            -15.10204082, -14.69387755, -14.28571429, -13.87755102,
            -13.46938776, -13.06122449, -12.65306122, -12.24489796,
            -11.83673469, -11.42857143, -11.02040816, -10.6122449 ,
            -10.20408163,  -9.79591837,  -9.3877551 ,  -8.97959184,
            -8.57142857,  -8.16326531,  -7.75510204,  -7.34693878,
            -6.93877551,  -6.53061224,  -6.12244898,  -5.71428571,
            -5.30612245,  -4.89795918,  -4.48979592,  -4.08163265,
            -3.67346939,  -3.26530612,  -2.85714286,  -2.44897959,
            -2.04081633,  -1.63265306,  -1.2244898 ,  -0.81632653,
            -0.40816327,   0.
    };
    for (float x: negFloatTests) {
        float stdExp = std::exp(x); // expected exponent
        float impExp = Exp(x);      // implemented exponent
        negFloatMaxErr = std::max(negFloatMaxErr, std::abs(stdExp - impExp));
    }

    float posFloatMaxErr = 0.0f; // relative error
    // Testing for x > 0
    float posFloatTests[] = {
            0.        ,  0.40816327,  0.81632653,  1.2244898 ,  1.63265306,
            2.04081633,  2.44897959,  2.85714286,  3.26530612,  3.67346939,
            4.08163265,  4.48979592,  4.89795918,  5.30612245,  5.71428571,
            6.12244898,  6.53061224,  6.93877551,  7.34693878,  7.75510204,
            8.16326531,  8.57142857,  8.97959184,  9.3877551 ,  9.79591837,
            10.20408163, 10.6122449 , 11.02040816, 11.42857143, 11.83673469,
            12.24489796, 12.65306122, 13.06122449, 13.46938776, 13.87755102,
            14.28571429, 14.69387755, 15.10204082, 15.51020408, 15.91836735,
            16.32653061, 16.73469388, 17.14285714, 17.55102041, 17.95918367,
            18.36734694, 18.7755102 , 19.18367347, 19.59183673, 20.
    };
    for (float x: posFloatTests) {
        float stdExp = std::exp(x); // expected exponent
        float impExp = Exp(x);      // implemented exponent
        float diff = std::abs(stdExp - impExp);
        posFloatMaxErr = std::max(posFloatMaxErr, diff / stdExp);
    }

    std::cout << "Max absolute error for x <= 0: " << negFloatMaxErr
    << "\nMax relative error for x > 0: " << posFloatMaxErr << '\n';
}
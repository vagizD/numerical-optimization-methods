#include "../include/diff.h"
#include "../include/stencil.h"

#include <cmath>
#include <iostream>
#include <functional>

double f(double x, double y) {
    return x * y + 5;
}

AAD22 g(AAD22 x, AAD22 y) {
    return cos(x) + y;
}


int main() {

    Differentiator<Derivative::X, DiffMethod::FwdADD, std::function<AAD22(AAD22, AAD22)>>(g, 5, 3);

    return 0;
}
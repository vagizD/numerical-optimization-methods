#pragma once
#include <iostream>

#include "container.h"

inline void test_container() {
    std::cout << "\n=================== START TEST CONTAINER ===================\n";

    std::cout << "TEST: METHODS & COPY CTOR\n";
    constexpr size_t N1 = 3;

    auto *A2D = new Matrix2D<N1>();
    A2D->init_zero();
    (*A2D)[5] = 1.5;

    Matrix2D<N1> *B2D = A2D;
    assert((*B2D)[5] == 1.5);

    Matrix1D<N1> A1D;
    A1D.init_zero();
    A1D[5] = -1.5;

    Matrix1D<N1> B1D = A1D;
    assert(B1D[5] == -1.5);

    constexpr size_t N2 = 100;

    MatrixSparse<N2, 5> ASparse;
    ASparse.init_zero();
    ASparse[10] = -1.5;

    MatrixSparse<N2, 5> BSparse = ASparse;
    assert(BSparse[10] == -1.5);

    std::cout << "PASSED\n\n";
    std::cout << "=================== END TEST CONTAINER ===================\n";
}


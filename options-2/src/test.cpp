#include "test_container.h"
#include "test_slae.h"
#include "test_stepper.h"
#include "test_utils.h"

int main() {
    test_container();
    test_utils();
    test_slae();
    test_stepper();
}

/*
           (V_i, S_j)       | Matrix2D | Matrix1D | MatrixGSL  |
Answer for (S_0, V_0) at t=T: 12.63876  | 12.63876  | 10.17400  |
32.10600 seconds

BSM:

*/

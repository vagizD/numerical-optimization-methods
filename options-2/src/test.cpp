#include "test_utils.h"
#include "test_container.h"
#include "test_slae.h"

int main() {
    test_container();
    test_utils();
    test_slae();
}

/* output

=================== START TEST SLAE ===================
TEST 1
x1 = 1, ans1 = 1
x2 = 2, ans2 = 2
x3 = 3, ans3 = 3
PASSED

TEST 2
x1 = 1, ans1 = 1
x2 = 0, ans2 = 0
x3 = 2, ans3 = 2
PASSED


==================== END TEST SLAE ====================


=================== START TEST UTILS ===================
TEST: ADJUST BOUNDARY
i          : 177
j          : 170
S_max      : 117.64706
S_0 * n / j: 117.64706
V_max      :   1.12994
V_0 * m / i:   1.12994
PASSED


==================== END TEST SLAE ====================

*/

#pragma once
#include <array>
#include <container.h>
#include <iostream>

// Given logical indexes (lrow, lcol), computes physical index (prow, pcol)
// and then returns physical offset. Here lrow == prow.
template <const size_t N>
size_t l2p(const size_t lrow, const size_t lcol, const std::array<size_t, N> &pofl) {
    return lrow * N + pofl[lcol];
}

// Find logical index of maximum absolute value in lrow-th row of matrix A,
// on interval [lcol_offset, N)
template <const size_t N, typename C>
size_t get_max_index(
    const C &A,
    const size_t lrow,
    const size_t lcol_offset,
    const std::array<size_t, N> &pofl
) {
    size_t ans = lcol_offset;
    for (size_t i = lcol_offset; i < N; ++i) {
        size_t ans_pindex = l2p<N>(lrow, ans, pofl);
        size_t cur_pindex = l2p<N>(lrow, i, pofl);
        if (std::abs(A[cur_pindex]) > std::abs(A[ans_pindex])) {
            ans = i;
        }
    }
    return ans;
}

// Solve SLAE Ax = b using Gaussian elimination
template <const size_t N, typename C>
void solve_slae_(
    C &A,
    std::array<double, N> &b,
    std::array<double, N> &x
) {
    // 1. Initialize supplementary mapping array
    std::array<size_t, N> pofl{};  // physical col number of logical col index
    for (size_t i = 0; i < N; ++i) {
        pofl[i] = i;
    }

    // 2. Transform A to Upper Triangular Matrix
    // Iterate by logical index of row/col
    for (size_t li = 0; li < N; ++li) {
        // Get logical index of maximum element in the row
        size_t lj = get_max_index(A, li, li, pofl);

        if (std::abs(A[l2p<N>(li, lj, pofl)]) < 1e-10) {
            std::cerr << "Zero as max element\n";
            exit(1);
        }

        // Swap logical indexes of leftmost and max columns
        std::swap(pofl[li], pofl[lj]);

        // Subtract current row from each row below
        size_t poffset = l2p<N>(li, li, pofl);  // A[i, i]
        for (size_t lrow = li + 1; lrow < N; ++lrow) {
            size_t row_poffset = l2p<N>(lrow, li, pofl);  //  A[i + r, i]
            if (A[row_poffset] == 0) {
                continue;
            }
            double m = A[row_poffset] / A[poffset];

            A[row_poffset] = 0;

            for (size_t lcol = li + 1; lcol < N; ++lcol) {
                size_t row_col_poffset = l2p<N>(lrow, lcol, pofl);  // A[i + r, i + c]
                size_t col_poffset = l2p<N>(li, lcol, pofl);        // A[i, i + c]
                A[row_col_poffset] -= m * A[col_poffset];
            }

            b[lrow] -= m * b[li];
        }
    }

    // 3. Backtrack and compute values of x

    // x_N = b_N / A[N, N];
    // x_i = (b_i - \sum_{j=i+1}^{N}x_j * A[i, j]) / A[i, i];
    // above is based on fact that x_j for j > i are already computed
    bool flag = true;
    for (size_t i = N - 1; flag; --i) {
        x[pofl[i]] = b[i];
        for (size_t j = i + 1; j < N; ++j) {
            x[pofl[i]] -= x[pofl[j]] * A[l2p<N>(i, j, pofl)];
        }
        x[pofl[i]] /= A[l2p<N>(i, i, pofl)];

        flag = (i != 0);
    }
}


// ================================ SPARSE MATRIX GE ================================
template <const size_t N>
size_t get_max_index_sparse(
    const MatrixSparse<N, 5> &A,
    const size_t lrow,
    const size_t lcol_offset,
    const std::array<size_t, N> &pofl
) {
    size_t ans = lcol_offset;
    for (size_t i = lcol_offset; i < std::min(lcol_offset + A.d + 1, N); ++i) {
        size_t ans_pindex = l2p<N>(lrow, ans, pofl);
        size_t cur_pindex = l2p<N>(lrow, i, pofl);
        if (std::abs(A[cur_pindex]) > std::abs(A[ans_pindex])) {
            ans = i;
        }
    }
    return ans;
}


template <const size_t N>
void solve_slae_sparse_(
    MatrixSparse<N, 5> &A,
    std::array<double, N> &b,
    std::array<double, N> &x
) {
    std::array<size_t, N> pofl{};  // physical col number of logical col index
    for (size_t i = 0; i < N; ++i) {
        pofl[i] = i;
    }

    for (size_t li = 0; li < N; ++li) {
        size_t lj = get_max_index_sparse(A, li, li, pofl);

        if (std::abs(A[l2p<N>(li, lj, pofl)]) < 1e-10) {
            std::cerr << "Zero as max element\n";
            exit(1);
        }

        std::swap(pofl[li], pofl[lj]);

        size_t poffset = l2p<N>(li, li, pofl);  // A[i, i]
        for (size_t lrow = li + 1; lrow < std::min(li + A.d + 1, N); ++lrow) {
            size_t row_poffset = l2p<N>(lrow, li, pofl);  //  A[i + r, i]
            double m = A[row_poffset] / A[poffset];

            A[row_poffset] = 0;

            for (size_t lcol = li + 1; lcol < std::min(li + A.d + 1, N); ++lcol) {
                size_t row_col_poffset = l2p<N>(lrow, lcol, pofl);  // A[i + r, i + c]
                size_t col_poffset = l2p<N>(li, lcol, pofl);        // A[i, i + c]
                A[row_col_poffset] -= m * A[col_poffset];
            }

            b[lrow] -= m * b[li];
        }
    }

    bool flag = true;
    for (size_t i = N - 1; flag; --i) {
        x[pofl[i]] = b[i];
        for (size_t j = i + 1; j < std::min(i + A.d + 1, N); ++j) {
            x[pofl[i]] -= x[pofl[j]] * A[l2p<N>(i, j, pofl)];
        }
        x[pofl[i]] /= A[l2p<N>(i, i, pofl)];

        flag = (i != 0);
    }
}

// ==================================================================================


template <const size_t N, typename C>
void solve_slae(
    const C &A,
    const std::array<double, N> &b,
    std::array<double, N> &x
) {
    std::array<double, N> b_copy = b;
    auto A_copy = A;
    solve_slae_(A_copy, b_copy, x);
}


template <const size_t N, typename C>
void solve_slae_sparse(
    const C &A,
    const std::array<double, N> &b,
    std::array<double, N> &x
) {
    std::array<double, N> b_copy = b;
    auto A_copy = A;
    solve_slae_sparse_(A_copy, b_copy, x);
}

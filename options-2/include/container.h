#pragma  once
#include <array>
#include <cassert>
#include <algorithm>

// map 2D index to 1D
inline size_t map2to1(const size_t i, const size_t j, const size_t n_cols) {
    return n_cols * i + j;
}

// Abstract class for any container for matrix A
class Container {
public:
    virtual double get(size_t i, size_t j) = 0;
    virtual void set(size_t i, size_t j, double x) = 0;
    virtual void add(size_t i, size_t j, double x) = 0;
    virtual void init_zero() = 0;
    virtual ~Container() = default;
};

// 2D matrix representation as array of arrays
// NOT a contiguous block of memory
template<const size_t N>
class Matrix2D final : Container {
    std::array<std::array<double, N>, N> A;

public:
    void init_zero() override {
        for (size_t i = 0; i < N; ++i) {
             std::fill(A[i].begin(), A[i].end(), 0);
        }
    }

    double get(const size_t i, const size_t j) override {
        return A[i][j];
    }

    void set(const size_t i, const size_t j, const double x) override {
        A[i][j] = x;
    }

    void add(const size_t i, const size_t j, const double x) override {
        A[i][j] += x;
    }
};


template<const size_t N>
class Matrix1D final : Container {
    std::array<double, N * N> A;

public:
    void init_zero() override {
        std::fill(A.begin(), A.end(), 0);
    }

    double get(const size_t i, const size_t j) override {
        return A[map2to1(i, j, N)];
    }

    void set(const size_t i, const size_t j, const double x) override {
        A[map2to1(i, j, N)] = x;
    }

    void add(const size_t i, const size_t j, const double x) override {
        A[map2to1(i, j, N)] += x;
    }
};


template<const size_t N>
class MatrixGSL final : Container {};


template<const size_t N>
class MatrixSparse final : Container {
    std::array<double, N * 5> A;

public:
    const size_t M = N * 5;

    // check that accessed elements exist in matrix
    static void check_boundary(const size_t i, const size_t j) {
        if (i >= j) {
            assert(i - j <= 2);
        } else {
            assert(j - i <= 2);
        }
    }

    void init_zero() override {
        std::fill(A.begin(), A.end(), 0);
    }

    double get(const size_t i, const size_t j) override {
        check_boundary(i, j);
        return A[map2to1(i, i - j + 2, 5)];
    }

    void set(const size_t i, const size_t j, const double x) override {
        check_boundary(i, j);
        A[map2to1(i, i - j + 2, 5)] = x;
    }

    void add(const size_t i, const size_t j, const double x) override {
        check_boundary(i, j);
        A[map2to1(i, i - j + 2, 5)] += x;
    }
};


// 11 12 13  0 ...
// 21 22 23 24  0 ...
// 31 32 33 34 35  0 ...
//  0 42 43 44 45 46  0 ...
//  0  0 53 54 55 56 57  0 ...
//  0  0  0 64 65 66 67 68  0 ...
//  0  0  0

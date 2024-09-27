#pragma  once
#include <array>
#include <cassert>
#include <algorithm>

enum class ContainerType {
    Matrix2D,
    Matrix1D,
    MatrixGSL,
    MatrixSparse
};

// map 2D index to 1D
inline size_t map2to1(const size_t i, const size_t j, const size_t n_cols) {
    return n_cols * i + j;
}

// map 1D index to 2D
inline std::pair<size_t, size_t> map1to2(const size_t offset, const size_t n_cols) {
    return {offset / n_cols, offset % n_cols};
}

// Abstract class for any container for matrix A
class Container {
public:
    Container() = default;

    [[nodiscard]] virtual ContainerType type() const = 0;
    [[nodiscard]] virtual const size_t size() const = 0;

    virtual double get(size_t i, size_t j) = 0;
    virtual void set(size_t i, size_t j, double x) = 0;
    virtual void add(size_t i, size_t j, double x) = 0;
    virtual void init_zero() = 0;
    virtual double operator[](size_t offset) const = 0;  // get method (if const ref)
    virtual double& operator[](size_t offset) = 0;       // set method
    virtual ~Container() = default;
};

// 2D matrix representation as array of arrays
// NOT a contiguous block of memory
template<const size_t N>
class Matrix2D final : Container {
    std::array<std::array<double, N>, N> A;

public:
    [[nodiscard]] ContainerType type() const override {
        return ContainerType::Matrix2D;
    }

    [[nodiscard]] const size_t size() const override {
        return N;
    }

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

    double operator[](const size_t offset) const override {
        auto [i, j] = map1to2(offset, size());
        return A[i][j];
    }

    double& operator[](const size_t offset) override {
        auto [i, j] = map1to2(offset, size());
        return A[i][j];
    }
};


template<const size_t N>
class Matrix1D final : Container {
    std::array<double, N * N> A;

public:
    [[nodiscard]] const size_t size() const override {
        return N;
    }

    [[nodiscard]] ContainerType type() const override {
        return ContainerType::Matrix1D;
    }

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

    double operator[](const size_t offset) const override {
        return A[offset];
    }

    double& operator[](const size_t offset) override {
        return A[offset];
    }
};


template<const size_t N>
class MatrixGSL final : Container {};


template<const size_t N, const size_t n_cols>
class MatrixSparse final : Container {
    std::array<double, N * n_cols> A;

public:
    const size_t d = n_cols / 2;
    const size_t M = N * n_cols;

    // check that accessed elements exist in matrix
    void check_boundary(const size_t i, const size_t j) const {
        if (i >= j) {
            assert(i - j <= d);
        } else {
            assert(j - i <= d);
        }
    }

    [[nodiscard]] const size_t size() const override {
        return N;
    }

    [[nodiscard]] ContainerType type() const override {
        return ContainerType::MatrixSparse;
    }

    void init_zero() override {
        std::fill(A.begin(), A.end(), 0);
    }

    double get(const size_t i, const size_t j) override {
        check_boundary(i, j);
        return A[map2to1(i, j + d - i, n_cols)];
    }

    void set(const size_t i, const size_t j, const double x) override {
        check_boundary(i, j);
        A[map2to1(i, j + d - i, n_cols)] = x;
    }

    void add(const size_t i, const size_t j, const double x) override {
        check_boundary(i, j);
        A[map2to1(i, j + d - i, n_cols)] += x;
    }

    double operator[](const size_t offset) const override {
        auto [i, j] = map1to2(offset, size());
        return A[map2to1(i, j + d - i, n_cols)];
    }

    double& operator[](const size_t offset) override {
        auto [i, j] = map1to2(offset, size());
        return A[map2to1(i, j + d - i, n_cols)];
    }
};


// 11 12 13  0 ...
// 21 22 23 24  0 ...
// 31 32 33 34 3n_cols  0 ...
//  0 42 43 44 4n_cols 46  0 ...
//  0  0 n_cols3 n_cols4 n_colsn_cols n_cols6 n_cols7  0 ...
//  0  0  0 64 6n_cols 66 67 68  0 ...
//  0  0  0

#pragma once

#include <cstddef>
#include <string>
#include <vector>

/// Floating-point precision
using Real = double;

/// Integer type
using Int  = int;

/// CGNS uses 64-bit integers for large grids
using CGNSSize = long long;

/// 3D vector
struct Vec3 {
    Real x, y, z;
    Vec3() : x(0), y(0), z(0) {}
    Vec3(Real x_, Real y_, Real z_) : x(x_), y(y_), z(z_) {}
};

/// 1D contiguous array with dynamically allocated storage
template <typename T>
class MultiArray1D {
public:
    MultiArray1D();
    MultiArray1D(std::size_t n);
    void allocate(std::size_t n);

    T&       operator[](std::size_t i);
    const T& operator[](std::size_t i) const;

    T*       data();
    const T* data() const;
    std::size_t size() const;
    void fill(T val);

private:
    std::size_t size_;
    std::vector<T> data_;
};

/// 3D contiguous array: index = i + ni*(j + nj*k)
template <typename T>
class MultiArray3D {
public:
    MultiArray3D();
    MultiArray3D(Int ni, Int nj, Int nk);
    void allocate(Int ni, Int nj, Int nk);

    T&       operator()(Int i, Int j, Int k);
    const T& operator()(Int i, Int j, Int k) const;

    T*       data();
    const T* data() const;
    std::size_t size() const;
    void fill(T val);

    Int ni() const;
    Int nj() const;
    Int nk() const;

private:
    Int ni_, nj_, nk_;
    std::vector<T> data_;
    std::size_t idx(Int i, Int j, Int k) const;
};

#include "types.hxx"

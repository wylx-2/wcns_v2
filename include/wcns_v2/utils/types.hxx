#pragma once

#include "types.h"
#include <algorithm>

// ---- MultiArray1D ----
template <typename T>
inline MultiArray1D<T>::MultiArray1D() : size_(0) {}

template <typename T>
inline MultiArray1D<T>::MultiArray1D(std::size_t n) : size_(n), data_(n) {}

template <typename T>
inline void MultiArray1D<T>::allocate(std::size_t n) { size_ = n; data_.resize(n); }

template <typename T>
inline T& MultiArray1D<T>::operator[](std::size_t i) { return data_[i]; }

template <typename T>
inline const T& MultiArray1D<T>::operator[](std::size_t i) const { return data_[i]; }

template <typename T>
inline T* MultiArray1D<T>::data() { return data_.data(); }

template <typename T>
inline const T* MultiArray1D<T>::data() const { return data_.data(); }

template <typename T>
inline std::size_t MultiArray1D<T>::size() const { return size_; }

template <typename T>
inline void MultiArray1D<T>::fill(T val) { std::fill(data_.begin(), data_.end(), val); }

// ---- MultiArray3D ----
template <typename T>
inline MultiArray3D<T>::MultiArray3D() : ni_(0), nj_(0), nk_(0) {}

template <typename T>
inline MultiArray3D<T>::MultiArray3D(Int ni, Int nj, Int nk)
    : ni_(ni), nj_(nj), nk_(nk), data_(static_cast<std::size_t>(ni) * nj * nk) {}

template <typename T>
inline void MultiArray3D<T>::allocate(Int ni, Int nj, Int nk) {
    ni_ = ni; nj_ = nj; nk_ = nk;
    data_.resize(static_cast<std::size_t>(ni) * nj * nk);
}

template <typename T>
inline T& MultiArray3D<T>::operator()(Int i, Int j, Int k) { return data_[idx(i,j,k)]; }

template <typename T>
inline const T& MultiArray3D<T>::operator()(Int i, Int j, Int k) const { return data_[idx(i,j,k)]; }

template <typename T>
inline T* MultiArray3D<T>::data() { return data_.data(); }

template <typename T>
inline const T* MultiArray3D<T>::data() const { return data_.data(); }

template <typename T>
inline std::size_t MultiArray3D<T>::size() const { return data_.size(); }

template <typename T>
inline void MultiArray3D<T>::fill(T val) { std::fill(data_.begin(), data_.end(), val); }

template <typename T>
inline Int MultiArray3D<T>::ni() const { return ni_; }
template <typename T>
inline Int MultiArray3D<T>::nj() const { return nj_; }
template <typename T>
inline Int MultiArray3D<T>::nk() const { return nk_; }

template <typename T>
inline std::size_t MultiArray3D<T>::idx(Int i, Int j, Int k) const {
    return static_cast<std::size_t>(i) + static_cast<std::size_t>(ni_) *
           (static_cast<std::size_t>(j) + static_cast<std::size_t>(nj_) *
            static_cast<std::size_t>(k));
}

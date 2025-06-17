//
// Created by iUV on 3/7/2025.
//
#include "geometry.h"

template <typename T>
 Matrix<T> operator* ( Matrix<T> const &cur,  Matrix4x4<T> const &other) {
    return cur.multiply(other);
}

template<typename T>
 Matrix<T> operator* ( Matrix4x4<T> const &cur,  Matrix<T> const &other) {
    return cur.multiply(other);
}

template class Matrix<float>;
template class Matrix4x4<float>;
template Matrix<float> operator*( Matrix<float> const & ,  Matrix4x4<float> const &);
template Matrix<float> operator*( Matrix4x4<float> const &,  Matrix<float> const &);
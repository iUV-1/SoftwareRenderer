//
// Created by iUV on 3/7/2025.
//
#include "geometry.h"



/// Create a homogonized matrix from a vector
Vec4f homogonize(Vec3f const &v, float h) {
    Vec4f result;
    result.x = v.x;
    result.y = v.y;
    result.z = v.z;
    result.w = h;
    return result;
}

/// De-homogonize it
Vec3f dehomogonize(Vec4f const &v) {
    if(v.w == 0.) {
        throw std::invalid_argument("what the fuck");
    }
    Vec3f result;
    result.x = v.x / v.w;
    result.y = v.y / v.w;
    result.z = v.z / v.w;
    return result;
}

///// Reflect a vector from a surface
////
//Vec3f reflect(Vec3f const n, Vec3f const i) {
//
//}

template<typename T>
Vec2<T>::Vec2(Matrix<T> mat) {
    if(mat.rows >= 2) {
        x = mat[0][0];
        y = mat[1][0];
        return;
    }
    if(mat.cols >= 2) {
        x = mat[0][0];
        y = mat[0][1];
        return;
    }
    throw std::invalid_argument("Matrix doesn't have the appropriate size to cast");
}

template<typename T>
Vec3<T>::Vec3(Matrix<T> mat) {
    if(mat.rows >= 3) {
        x = mat[0][0];
        y = mat[1][0];
        z = mat[2][0];
        return;
    }
    if(mat.cols >= 3) {
        x = mat[0][0];
        y = mat[0][1];
        z = mat[0][2];
        return;
    }
    throw std::invalid_argument("Matrix doesn't have the appropriate size to cast");
}

template<typename T>
Vec4<T>::Vec4(Matrix<T> mat) {
    if(mat.rows >= 4) {
        x = mat[0][0];
        y = mat[1][0];
        z = mat[2][0];
        w = mat[3][0];
        return;
    }
    if(mat.cols >= 4) {
        x = mat[0][0];
        y = mat[0][1];
        z = mat[0][2];
        w = mat[0][3];
        return;
    }
    throw std::invalid_argument("Matrix doesn't have the appropriate size to cast");
}

template<typename T>
Vec4<T> Matrix4x4<T>::operator*(const Vec4<T> &other) {
    Vec4<T> result;
    const auto& A = this->data;
    for(int i = 0; i < 4; i++) {
        result[i] = other.x * A[i][0] + other.y * A[i][1] + other.z * A[i][2] + other.w * A[i][3];
    }
    return result;
}

template<typename T>
Vec3<T> Matrix3x3<T>::operator*(const Vec3<T> &other) {
    Vec3<T> result;
    const auto& A = this->data;
#pragma unroll
    for(int i = 0; i < 3; i++) {
        result[i] = other.x * A[i][0] + other.y * A[i][1] + other.z * A[i][2];
    }
    return result;
}

template<typename T>
Matrix<T> Matrix<T>::operator*(const Vec2<T> &other) {
    if(cols != 2) {
        throw std::invalid_argument("Invalid matrix size");
    }

    Matrix<T> result(rows, 1);
    const auto& A = this->data;
    for(int i = 0; i < rows; i++) {
        result[i][0] = other.x * A[i][0] + other.y * A[i][1];
    }
    return result;
}

template<typename T>
Matrix<T> Matrix<T>::operator*(const Vec3<T> &other) {
    if(cols != 3) {
        throw std::invalid_argument("Invalid matrix size");
    }
    Matrix<T> result(rows, 1);
    const auto& A = this->data;
    for(int i = 0; i < rows; i++) {
        result[i][0] = other.x * A[i][0] + other.y * A[i][1] + other.z * A[i][2];
    }
    return result;
}
template class Vec2<float>;
template class Vec3<float>;
template class Matrix<float>;
template class Matrix4x4<float>;
template class Matrix3x3<float>;
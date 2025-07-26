 //
// Created by iUV on 9/7/2024.™
//

#ifndef SOFTWARERENDERER_GEOMETRY_H
#define SOFTWARERENDERER_GEOMETRY_H

#endif //SOFTWARERENDERER_GEOMETRY_H
#ifndef __GEOMETRY_H__
#define __GEOMETRY_H__

#include <cmath>
#include <vector>
#include <iostream>

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

template <class t> struct Vec2 {
    union {
        struct {t u, v;};
        struct {t x, y;};
        t raw[2];
    };
    Vec2() : u(0), v(0) {}
    Vec2(t _u, t _v) : u(_u),v(_v) {}
    inline Vec2<t> operator +(const Vec2<t> &V) const { return Vec2<t>(u+V.u, v+V.v); }
    inline Vec2<t> operator -(const Vec2<t> &V) const { return Vec2<t>(u-V.u, v-V.v); }
    inline Vec2<t> operator *(float f)          const { return Vec2<t>(u*f, v*f); }
    template <class > friend std::ostream& operator<<(std::ostream& s, Vec2<t>& v);
    t& operator[](size_t const &index) {
        if(index == 0) return x;
        if(index == 1) return y;
        throw std::out_of_range("index out of range");
    };
};

template <class t> struct Vec3 {
    union {
        struct {t x, y, z;};
        struct { t ivert, iuv, inorm; };
        t raw[3];
    };
    Vec3() : x(0), y(0), z(0) {}
    Vec3(t _x, t _y, t _z) : x(_x),y(_y),z(_z) {}

    inline Vec3<t> operator ^(const Vec3<t> &v) const { return Vec3<t>(y*v.z-z*v.y, z*v.x-x*v.z, x*v.y-y*v.x); }
    inline Vec3<t> operator +(const Vec3<t> &v) const { return Vec3<t>(x+v.x, y+v.y, z+v.z); }
    inline Vec3<t> operator -(const Vec3<t> &v) const { return Vec3<t>(x-v.x, y-v.y, z-v.z); }
    inline Vec3<t> operator *(float f)          const { return Vec3<t>(x*f, y*f, z*f); }
    inline t       operator *(const Vec3<t> &v) const { return x*v.x + y*v.y + z*v.z; }
    float norm () const { return std::sqrt(x*x+y*y+z*z); }
    Vec3<t> & normalize(t l=1) { *this = (*this)*(l/norm()); return *this; }
    template <class > friend std::ostream& operator<<(std::ostream& s, Vec3<t>& v);
    t& operator[](size_t const &index) {
        if(index == 0) return x;
        if(index == 1) return y;
        if(index == 2) return z;
        throw std::out_of_range("index out of range");
    };
};

typedef Vec2<float> Vec2f;
typedef Vec2<int>   Vec2i;
typedef Vec3<float> Vec3f;
typedef Vec3<int>   Vec3i;

template <class t> std::ostream& operator<<(std::ostream& s, Vec2<t>& v) {
    s << "(" << v.x << ", " << v.y << ")\n";
    return s;
}

template <class t> std::ostream& operator<<(std::ostream& s, Vec3<t>& v) {
    s << "(" << v.x << ", " << v.y << ", " << v.z << ")\n";
    return s;
}

template <typename T>
class Matrix4x4;

template <typename T> class Matrix {
public:
    std::vector<std::vector<T>> data;
    size_t rows, cols;

    Matrix(size_t rows, size_t columns, T defaultVal = T{}): rows(rows), cols(columns),data(rows, std::vector<T>(columns, defaultVal)) {}

    // init a row Matrix from a 3D vector
    explicit Matrix(Vec3<T> vec): Matrix(3, 1){
        data[0][0] = vec.x;
        data[1][0] = vec.y;
        data[2][0] = vec.z;
    }
    std::vector<T>& operator[](size_t row) {return data[row];}
    const std::vector<T>& operator[](size_t row) const {return data[row];} // return a const reference to a row

    Matrix<T> multiply(Matrix<T> matrix) const {
        if (cols != matrix.rows) {
            std::cerr << "Matrix multiplication mismatch!!" << std::endl;
            // unable to do it, what should i return?
            return matrix;
        }

        Matrix<T> result = Matrix(rows, matrix.cols);
        for(int i = 0; i < result.rows; i++) {
            for (int j = 0; j < result.cols; j++){
                for (int k = 0; k < cols; k++) {
                    result[i][j] += data[i][k] * matrix[k][j];
                }
            }
        }
        return result;
    }

    void display() const {
        for (size_t i = 0; i < rows; i++) {
            for (size_t j = 0; j < cols; j++) {
                std::cout << data[i][j] << " ";
            }
            std::cout << std::endl;
        }
    }



    Matrix<T> operator*(Matrix<T> const &other) { return this->multiply(other); }
    template<class U>
    friend Matrix<U> operator*( Matrix const &cur,  Matrix4x4<U> const &other);
};

template <class T> std::ostream& operator<<(std::ostream& os, Matrix<T> m) {
    for (size_t i = 0; i < m.rows; i++) {
        for (size_t j = 0; j < m.cols; j++) {
            os << m[i][j] << " ";
        }
        os << "\n";
    }
    return os;
}

// 3x3 matrix derived from Matrix class
template <typename T> class Matrix3x3: public Matrix<T> {
public:
    Matrix3x3() : Matrix<T>(3, 3) {}

    // optimized calculation just for 3x3
    Matrix3x3 multiply3x3(const Matrix3x3 &other) {
        Matrix3x3 result;

        const auto& A = this->data;
        const auto& B = other.data;

        result(0, 0) = A[0][0] * B[0][0] + A[0][1] * B[1][0] + A[0][2] * B[2][0];
        result(0, 1) = A[0][0] * B[0][1] + A[0][1] * B[1][1] + A[0][2] * B[2][1];
        result(0, 2) = A[0][0] * B[0][2] + A[0][1] * B[1][2] + A[0][2] * B[2][2];

        result(1, 0) = A[1][0] * B[0][0] + A[1][1] * B[1][0] + A[1][2] * B[2][0];
        result(1, 1) = A[1][0] * B[0][1] + A[1][1] * B[1][1] + A[1][2] * B[2][1];
        result(1, 2) = A[1][0] * B[0][2] + A[1][1] * B[1][2] + A[1][2] * B[2][2];

        result(2, 0) = A[2][0] * B[0][0] + A[2][1] * B[1][0] + A[2][2] * B[2][0];
        result(2, 1) = A[2][0] * B[0][1] + A[2][1] * B[1][1] + A[2][2] * B[2][1];
        result(2, 2) = A[2][0] * B[0][2] + A[2][1] * B[1][2] + A[2][2] * B[2][2];

        return result;
    }
};

typedef Matrix3x3<float> Matrix3x3if;

// 4x4 matrix derived from Matrix class
template <typename T> class Matrix4x4: public Matrix<T> {
public:
     Matrix4x4() : Matrix<T>(4, 4) {}

     // optimized calculation just for 4x4
     Matrix4x4 multiply4x4(const Matrix4x4 &other) {
         Matrix4x4 result;

         const auto& A = this->data;
         const auto& B = other.data;

         for (int i = 0; i < 4; i++) {
             result[i][0] = A[i][0] * B[0][0] + A[i][1] * B[1][0] + A[i][2] * B[2][0] + A[i][3] * B[3][0];
             result[i][1] = A[i][0] * B[0][1] + A[i][1] * B[1][1] + A[i][2] * B[2][1] + A[i][3] * B[3][1];
             result[i][2] = A[i][0] * B[0][2] + A[i][1] * B[1][2] + A[i][2] * B[2][2] + A[i][3] * B[3][2];
             result[i][3] = A[i][0] * B[0][3] + A[i][1] * B[1][3] + A[i][2] * B[2][3] + A[i][3] * B[3][3];
         }

         return result;
     }

     // create a 4x4 identity matrix
     static Matrix4x4 identity() {
         Matrix4x4 result;

         for (int i = 4; i--; ) {
             result[i][i] = 1;
         }
         return result;
     }

     // Invert its own matrix
//     void invert() {
//        Matrix4x4<T> temp = *this;
//
//        T det = 0;
//
//        // diagonal
//        T a = 1;
//        T b = 1;
//        for(int i = 0; i < 4; ++i) {
//            a *= this[i][i];
//            b *= this[4 - i][4 - i];
//        }
//     }

    // From ChatGPT lol
    void inverseTranspose() {
        Matrix4x4<T> temp = *this;

        // Step 1: Compute determinant via cofactor expansion along row 0
        T det = 0;
        for (int j = 0; j < 4; ++j) {
            Matrix<T> minor(3, 3);
            for (int mi = 1; mi < 4; ++mi) {
                int r = mi - 1;
                int c = 0;
                for (int mj = 0; mj < 4; ++mj) {
                    if (mj == j) continue;
                    minor[r][c++] = temp[mi][mj];
                }
            }

            T cofactor =
                    minor[0][0]*(minor[1][1]*minor[2][2] - minor[1][2]*minor[2][1]) -
                    minor[0][1]*(minor[1][0]*minor[2][2] - minor[1][2]*minor[2][0]) +
                    minor[0][2]*(minor[1][0]*minor[2][1] - minor[1][1]*minor[2][0]);

            det += ((j % 2 == 0 ? 1 : -1) * temp[0][j] * cofactor);
        }

        if (std::abs(det) < std::numeric_limits<T>::epsilon() * std::abs(det)) {
            std::cerr << "Matrix is not invertible!\n";
            *this = Matrix4x4<T>::identity();
            return;
        }

        // Step 2: Compute inverse and transpose at the same time (adjugate / det)
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                Matrix<T> minor(3, 3);
                for (int mi = 0, r = 0; mi < 4; mi++) {
                    if (mi == i) continue;
                    for (int mj = 0, c = 0; mj < 4; mj++) {
                        if (mj == j) continue;
                        minor[r][c++] = temp[mi][mj];
                    }
                    r++;
                }

                T cofactor =
                        minor[0][0]*(minor[1][1]*minor[2][2] - minor[1][2]*minor[2][1]) -
                        minor[0][1]*(minor[1][0]*minor[2][2] - minor[1][2]*minor[2][0]) +
                        minor[0][2]*(minor[1][0]*minor[2][1] - minor[1][1]*minor[2][0]);

                (*this)[j][i] = ((i + j) % 2 == 0 ? 1 : -1) * cofactor / det;  // transpose here
            }
        }
    }


     inline Matrix4x4<T> operator*( Matrix4x4<T> const &other) { return this->multiply4x4(other); };
     template <typename U>
     friend  Matrix<U> operator* ( Matrix4x4<U> const &cur,  Matrix<U> const &other);
};

 template <typename T>
 Matrix<T> operator* ( Matrix<T> const &cur,  Matrix4x4<T> const &other) {
     return cur.multiply(other);
 }

 template<typename T>
 Matrix<T> operator* ( Matrix4x4<T> const &cur,  Matrix<T> const &other) {
     return cur.multiply(other);
 }
 template Matrix<float> operator*( Matrix<float> const & ,  Matrix4x4<float> const &);
 template Matrix<float> operator*( Matrix4x4<float> const &,  Matrix<float> const &);
typedef Matrix4x4<float> Matrix4x4f;
#endif //__GEOMETRY_H__
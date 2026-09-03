 //
// Created by iUV on 9/7/2024.™
//

#pragma once
#include <cmath>
#include <vector>
#include <iostream>
#include <array>

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

constexpr float MAX_FLOAT = std::numeric_limits<float>::max();
constexpr float EPS = std::numeric_limits<float>::epsilon();

// Forward declaration for Vectors
template <typename T, size_t w, size_t h>
class Matrix;

template <class t> class Vec2 {
public:
    union {
        struct {t u, v;};
        struct {t x, y;};
        t raw[2];
    };
    Vec2() : u(0), v(0) {}
    Vec2(t _u, t _v) : u(_u),v(_v) {}
    Vec2(Matrix<t, 2, 1> mat) {
        x = mat(0,0);
        y = mat(0,1);
    }
    Vec2(Matrix<t, 1, 2> mat) {
        x = mat(0, 0);
        y = mat(1, 0);
    }

//    Vec2(Matrix<t, > mat);
    inline Vec2<t> operator +(const Vec2<t> &V) const { return Vec2<t>(u+V.u, v+V.v); }
    inline Vec2<t> operator -(const Vec2<t> &V) const { return Vec2<t>(u-V.u, v-V.v); }
    inline Vec2<t> operator *(float f)          const { return Vec2<t>(u*f, v*f); }
    template <class > friend std::ostream& operator<<(std::ostream& s, Vec2<t>& v);
    Vec2<t> & normalize(t l=1) { *this = (*this)*(l/norm()); return *this; }
    float norm () const { return std::sqrt(x*x+y*y); }
    t& operator[](size_t const &index) {
        if(index == 0) return x;
        if(index == 1) return y;
        throw std::out_of_range("index out of range");
    };
};

// Forward declaration
// For homogonization implementation
template<class u>
class Vec4;

template <class t> class Vec3 {
public:
    union {
        struct {t x, y, z;};
        struct { t ivert, iuv, inorm; };
        t raw[3];
    };
    Vec3() : x(0), y(0), z(0) {}
    Vec3(t _x, t _y, t _z) : x(_x),y(_y),z(_z) {}
    Vec3(Matrix<t, 3, 1> mat) {
        x = mat(0, 0);
        y = mat(0, 1);
        z = mat(0, 2);
    }

    Vec3(Matrix<t, 1, 3> mat) {
        x = mat(0, 0);
        y = mat(1, 0);
        z = mat(2, 0);
    }
    // Dehomogonize
    Vec3(Vec4<t> other): x(other.x / (other.w + EPS)), y(other.y / (other.w + EPS)), z(other.z / (other.w + EPS)) {}
    ~Vec3() {}

    inline Vec3<t> operator ^(const Vec3<t> &v) const { return Vec3<t>(y*v.z-z*v.y, z*v.x-x*v.z, x*v.y-y*v.x); }
    inline Vec3<t> operator +(const Vec3<t> &v) const { return Vec3<t>(x+v.x, y+v.y, z+v.z); }
    inline Vec3<t> operator -(const Vec3<t> &v) const { return Vec3<t>(x-v.x, y-v.y, z-v.z); }
    inline Vec3<t> operator *(float f)          const { return Vec3<t>(x*f, y*f, z*f); }
    inline t       operator *(const Vec3<t> &v) const { return x*v.x + y*v.y + z*v.z; }
    float norm () const { return std::sqrt(x*x+y*y+z*z); }
    Vec3<t> & normalize(t l=1) {

        *this = (*this) * (l/norm());
        return *this;
    }
    template <class > friend std::ostream& operator<<(std::ostream& s, Vec3<t>& v);
    t& operator[](size_t const &index) {
        if(index == 0) return x;
        if(index == 1) return y;
        if(index == 2) return z;
        throw std::out_of_range("index out of range");
    }

    Vec3<float> &operator*=(float a) {
        x *= a;
        y *= a;
        z *= a;
        return *this;
    };

    // Move constructor
    Vec3(Vec3&& other) noexcept {
        x = other.x;
        y = other.y;
        z = other.z;
    }
    // Copy constructor
    Vec3(const Vec3 &other) {
        x = other.x;
        y = other.y;
        z = other.z;
    }
    // Copy operator
    Vec3& operator=(const Vec3 &other) {
        x = other.x;
        y = other.y;
        z = other.z;
        return *this;
    }
};

template <class t> class Vec4 {
public:
    union {
     struct {t x, y, z, w;};
     struct { t ivert, iuv, inorm, isomething; };
     t raw[4];
    };
    Vec4() : x(0), y(0), z(0), w(0) {}
    Vec4(t _x, t _y, t _z, t _w) : x(_x),y(_y),z(_z),w(_z) {}
    Vec4(Matrix<t, 4, 1> mat): x(mat(0, 0)), y(mat(1, 0)), z(mat(2, 0)), w(mat(3, 0)) { }
    Vec4(Matrix<t, 1, 4> mat): x(mat(0, 0)), y(mat(0, 1)), z(mat(0, 2)), w(mat(0, 3)) { }

    Vec4(Vec3<t> v, t h): x(v.x), y(v.y), z(v.z), w(h) {}

    inline Vec4<t> operator ^(const Vec4<t> &v) const { return Vec4<t>(y*v.z-z*v.y, z*v.x-x*v.z, x*v.y-y*v.x); }
    inline Vec4<t> operator +(const Vec4<t> &v) const { return Vec4<t>(x+v.x, y+v.y, z+v.z); }
    inline Vec4<t> operator -(const Vec4<t> &v) const { return Vec4<t>(x-v.x, y-v.y, z-v.z); }
    inline Vec4<t> operator *(float f)          const { return Vec4<t>(x*f, y*f, z*f); }
    inline t       operator *(const Vec4<t> &v) const { return x*v.x + y*v.y + z*v.z + w*v.w; }
    float norm () const { return std::sqrt(x*x+y*y+z*z+w*w); }
    Vec4<t> & normalize(t l=1) { *this = (*this)*(l/norm()); return *this; }
    template <class > friend std::ostream& operator<<(std::ostream& s, Vec4<t>& v);
    t& operator[](size_t const &index) {
     if(index == 0) return x;
     if(index == 1) return y;
     if(index == 2) return z;
     if(index == 3) return w;
     throw std::out_of_range("index out of range");
    }

    Vec4<t> &operator*=(t a) {
     x *= a;
     y *= a;
     z *= a;
     w *= a;
     return *this;
    };
};

typedef Vec2<float> Vec2f;
typedef Vec2<int>   Vec2i;
typedef Vec3<float> Vec3f;
typedef Vec3<int>   Vec3i;
typedef Vec4<int>   Vec4i;
typedef Vec4<float> Vec4f;

template <class t> std::ostream& operator<<(std::ostream& s, Vec2<t>& v) {
    s << "(" << v.x << ", " << v.y << ")\n";
    return s;
}

template <class t> std::ostream& operator<<(std::ostream& s, Vec3<t>& v) {
    s << "(" << v.x << ", " << v.y << ", " << v.z << ")\n";
    return s;
}

template <class t> std::ostream& operator<<(std::ostream& s, Vec4<t>& v) {
    s << "(" << v.x << ", " << v.y << ", " << v.z << ", " << v.w << ")\n";
    return s;
}

template <typename T>
class Matrix4x4;

// Matrix is row major
template <typename T, size_t w, size_t h>
class Matrix {
protected:
    /* Data */
    size_t rows, cols;
    std::array<T, w*h> data {};
public:
    /* Initializers */
    constexpr Matrix(): rows(w), cols(h) {}
    constexpr Matrix(std::initializer_list<T> list): rows(w), cols(h) {
        std::copy(list.begin(), list.end(), data.begin());
    }
    // init a row Matrix from a 3D vector
    explicit Matrix(Vec3<T> vec) {
        data[0] = vec.x;
        data[1 * rows + 0] = vec.y;
        data[2 * rows + 0] = vec.z;
    }

    /* Operator overloads */
    template<size_t other_w>
    Matrix<T, other_w, h> operator*(Matrix<T, other_w, w> const &other) { return this->multiply(other); }

    T &operator() (size_t x, size_t y) {
#ifndef NDEBUG
        if(x >= rows) {
            std::cerr << "Out of range access for x\n";
        }

        if(y >= cols) {
            std::cerr << "Out of range access for y\n";
        }
#endif
        return data[y * rows + x];
    }
    const T &operator() (size_t x, size_t y) const { return data[y * rows + x]; }

    Matrix<T, w, 1> operator*(Vec2<T> const &other) {
#ifndef NDEBUG
        if(cols != 2) {
            throw std::invalid_argument("Invalid matrix size");
        }
#endif
        Matrix<T, w, 1> result;
        for(int i = 0; i < rows; i++) {
            result(0, i) = other.x * data[0 * rows + i] + other.y * data[1 * rows + i];
        }
        return result;
    }

    Matrix<T, w, 1> operator*(Vec3<T> const &other) {
#ifndef NDEBUG
        if(cols != 3) {
            throw std::invalid_argument("Invalid matrix size");
        }
#endif
        Matrix<T, w, 1> result;
        for(int i = 0; i < rows; i++) {
            result(0, i) = other.x * data[0 * rows + 1] + other.y * data[1 * rows + i] + other.z * data[2 * rows + i];
        }
        return result;
    }

    // Not using it, couldn't be bothered to implement.
    // Matrix<T, w, h> operator*(Vec4<T> const &other);

    /* Functions */
    template<size_t other_w>
    Matrix<T, other_w, h> multiply(Matrix<T, other_w, w> matrix) const {
        Matrix<T, other_w, h> result;
        for(int i = 0; i < other_w; i++) {
            for (int j = 0; j < h; j++){
                for (int k = 0; k < h; k++) {
                    result(i, j) += data[i * rows + k] * matrix(k, j);
                }
            }
        }
        return result;
    }

    void display() const {
        for (size_t i = 0; i < rows; i++) {
            for (size_t j = 0; j < cols; j++) {
                std::cout << data[i * rows + j] << " ";
            }
            std::cout << "\n";
        }
    }

    void set_col(int idx, std::vector<T> vec) {
#ifndef NDEBUG
        if(vec.size() > cols) {
            std::cerr << "There are more elements in the vector than available columns!" << "\n";
        }
#endif
        for(int i = 0; i < cols; i++) {
            data[i][idx] = vec[i];
        }
    }

    void set_col(int idx, Vec2<T> v) {
        data[0 * rows + idx] = v.x;
        data[1 * rows + idx] = v.y;
    }

    void set_col(int idx, Vec3<T> v) {
        data[0 * rows + idx] = v.x;
        data[1 * rows + idx] = v.y;
        data[2 * rows + idx] = v.z;
    }

    void set_row(int idx, std::vector<T> vec) {
        if(vec.size() > rows) {
            std::cerr << "There are more elements in the vector than available rows!" << "\n";
        }
        for(int i = 0; i < cols; i++) {
            data[idx][i] = vec[i];
        }
    }

    void set_row(int idx, Vec2<T> v) {
        data[idx][0] = v.x;
        data[idx][1] = v.y;
    }

    void set_row(int idx, Vec3<T> v) {
        data[idx][0] = v.x;
        data[idx][1] = v.y;
        data[idx][2] = v.z;
    }
};

template <class T, size_t w, size_t h>
std::ostream& operator<<(std::ostream& os, Matrix<T, w ,h> m) {
    for (size_t i = 0; i < m.rows; i++) {
        for (size_t j = 0; j < m.cols; j++) {
            os << m(i,j) << " ";
        }
        os << "\n";
    }
    return os;
}

// 3x3 matrix derived from Matrix class
template <typename T>
class Matrix3x3: public Matrix<T, 3, 3> {
public:
    Matrix3x3() : Matrix<T, 3, 3>() {}

    /* Operator overloads */
    Vec3<T> operator*(const Vec3<T> &other) {
        Vec3<T> result;
        for(int i = 0; i < 3; i++) {
            result[i] = other.x * this->data[i * this->rows + 0]
                        + other.y * this->data[i * this->rows + 1]
                        + other.z * this->data[i * this->rows + 2];
        }
        return result;
    }

    // optimized calculation just for 3x3
//    Matrix3x3 multiply3x3(const Matrix3x3 &other) {
//        Matrix3x3 result;
//
//        const auto& A = this->data;
//        const auto& B = other.data;
//
//        for (int i = 0; i < 3; i++) {
//            result[i][0] = A[i][0] * B[0][0] + A[i][1] * B[1][0] + A[i][2] * B[2][0];
//            result[i][1] = A[i][0] * B[0][1] + A[i][1] * B[1][1] + A[i][2] * B[2][1];
//            result[i][2] = A[i][0] * B[0][2] + A[i][1] * B[1][2] + A[i][2] * B[2][2];
//        }
//
//        return result;
//    }
};

typedef Matrix3x3<float> Matrix3x3f;

// 4x4 matrix derived from Matrix class
template <typename T>
class Matrix4x4: public Matrix<T, 4, 4> {
public:
    constexpr Matrix4x4() {}
    constexpr Matrix4x4(std::initializer_list<T> list): Matrix<T, 4, 4>(list) {}
    // optimized calculation just for 4x4
    Matrix4x4 multiply4x4(const Matrix4x4 &other) {
        Matrix4x4 result;

        for (int i = 0; i < 4; i++) {
            // Use (*this)(row, col) to access the current matrix
            result(0, i) = (*this)(0, i) * other(0, 0) + (*this)(1, i) * other(0, 1) + (*this)(2, i) * other(0, 2) + (*this)(3, i) * other(0, 3);
            result(1, i) = (*this)(0, i) * other(1, 0) + (*this)(1, i) * other(1, 1) + (*this)(2, i) * other(1, 2) + (*this)(3, i) * other(1, 3);
            result(2, i) = (*this)(0, i) * other(2, 0) + (*this)(1, i) * other(2, 1) + (*this)(2, i) * other(2, 2) + (*this)(3, i) * other(2, 3);
            result(3, i) = (*this)(0, i) * other(3, 0) + (*this)(1, i) * other(3, 1) + (*this)(2, i) * other(3, 2) + (*this)(3, i) * other(3, 3);
        }

        return result;
    }

    // create a 4x4 identity matrix
    constexpr static Matrix4x4 create_identity() {
        Matrix4x4 result;

        for (int i = 4; i--; ) {
            result(i, i) = 1;
        }
        return result;
    }

    // Invert its own matrix
    // from chatGPT
    void invert() {
        Matrix4x4<T> inv;
        const auto &m = (*this);

        inv(0, 0) =  m(0, 1)*m(1, 2)*m(2, 3) - m(0, 1)*m(1, 3)*m(2, 2) - m(1, 1)*m(0, 2)*m(2, 3)
                     + m(1, 1)*m(0, 3)*m(2, 2) + m(2, 1)*m(0, 2)*m(1, 3) - m(2, 1)*m(0, 3)*m(1, 2);
        inv(1, 0) = -m(1, 0)*m(2, 2)*m(3, 3) + m(1, 0)*m(2, 3)*m(3, 2) + m(2, 0)*m(1, 2)*m(3, 3)
                    - m(2, 0)*m(1, 3)*m(3, 2) - m(3, 0)*m(1, 2)*m(2, 3) + m(3, 0)*m(1, 3)*m(2, 2);
        inv(2, 0) =  m(1, 0)*m(2, 1)*m(3, 3) - m(1, 0)*m(2, 3)*m(3, 1) - m(2, 0)*m(1, 1)*m(3, 3)
                     + m(2, 0)*m(1, 3)*m(3, 1) + m(3, 0)*m(1, 1)*m(2, 3) - m(3, 0)*m(1, 3)*m(2, 1);
        inv(3, 0) = -m(1, 0)*m(2, 1)*m(3, 2) + m(1, 0)*m(2, 2)*m(3, 1) + m(2, 0)*m(1, 1)*m(3, 2)
                    - m(2, 0)*m(1, 2)*m(3, 1) - m(3, 0)*m(1, 1)*m(2, 2) + m(3, 0)*m(1, 2)*m(2, 1);

        inv(0, 1) = -m(0, 1)*m(2, 2)*m(3, 3) + m(0, 1)*m(2, 3)*m(3, 2) + m(2, 1)*m(0, 2)*m(3, 3)
                    - m(2, 1)*m(0, 3)*m(3, 2) - m(3, 1)*m(0, 2)*m(2, 3) + m(3, 1)*m(0, 3)*m(2, 2);
        inv(1, 1) =  m(0, 0)*m(2, 2)*m(3, 3) - m(0, 0)*m(2, 3)*m(3, 2) - m(2, 0)*m(0, 2)*m(3, 3)
                     + m(2, 0)*m(0, 3)*m(3, 2) + m(3, 0)*m(0, 2)*m(2, 3) - m(3, 0)*m(0, 3)*m(2, 2);
        inv(2, 1) = -m(0, 0)*m(2, 1)*m(3, 3) + m(0, 0)*m(2, 3)*m(3, 1) + m(2, 0)*m(0, 1)*m(3, 3)
                    - m(2, 0)*m(0, 3)*m(3, 1) - m(3, 0)*m(0, 1)*m(2, 3) + m(3, 0)*m(0, 3)*m(2, 1);
        inv(3, 1) =  m(0, 0)*m(2, 1)*m(3, 2) - m(0, 0)*m(2, 2)*m(3, 1) - m(2, 0)*m(0, 1)*m(3, 2)
                     + m(2, 0)*m(0, 2)*m(3, 1) + m(3, 0)*m(0, 1)*m(2, 2) - m(3, 0)*m(0, 2)*m(2, 1);

        inv(0, 2) =  m(0, 1)*m(1, 2)*m(3, 3) - m(0, 1)*m(1, 3)*m(3, 2) - m(1, 1)*m(0, 2)*m(3, 3)
                     + m(1, 1)*m(0, 3)*m(3, 2) + m(3, 1)*m(0, 2)*m(1, 3) - m(3, 1)*m(0, 3)*m(1, 2);
        inv(1, 2) = -m(0, 0)*m(1, 2)*m(3, 3) + m(0, 0)*m(1, 3)*m(3, 2) + m(1, 0)*m(0, 2)*m(3, 3)
                    - m(1, 0)*m(0, 3)*m(3, 2) - m(3, 0)*m(0, 2)*m(1, 3) + m(3, 0)*m(0, 3)*m(1, 2);
        inv(2, 2) =  m(0, 0)*m(1, 1)*m(3, 3) - m(0, 0)*m(1, 3)*m(3, 1) - m(1, 0)*m(0, 1)*m(3, 3)
                     + m(1, 0)*m(0, 3)*m(3, 1) + m(3, 0)*m(0, 1)*m(1, 3) - m(3, 0)*m(0, 3)*m(1, 1);
        inv(3, 2) = -m(0, 0)*m(1, 1)*m(3, 2) + m(0, 0)*m(1, 2)*m(3, 1) + m(1, 0)*m(0, 1)*m(3, 2)
                    - m(1, 0)*m(0, 2)*m(3, 1) - m(3, 0)*m(0, 1)*m(1, 2) + m(3, 0)*m(0, 2)*m(1, 1);

        inv(0, 3) = -m(0, 1)*m(1, 2)*m(2, 3) + m(0, 1)*m(1, 3)*m(2, 2) + m(1, 1)*m(0, 2)*m(2, 3)
                    - m(1, 1)*m(0, 3)*m(2, 2) - m(2, 1)*m(0, 2)*m(1, 3) + m(2, 1)*m(0, 3)*m(1, 2);
        inv(1, 3) =  m(0, 0)*m(1, 2)*m(2, 3) - m(0, 0)*m(1, 3)*m(2, 2) - m(1, 0)*m(0, 2)*m(2, 3)
                     + m(1, 0)*m(0, 3)*m(2, 2) + m(2, 0)*m(0, 2)*m(1, 3) - m(2, 0)*m(0, 3)*m(1, 2);
        inv(2, 3) = -m(0, 0)*m(1, 1)*m(2, 3) + m(0, 0)*m(1, 3)*m(2, 1) + m(1, 0)*m(0, 1)*m(2, 3)
                    - m(1, 0)*m(0, 3)*m(2, 1) - m(2, 0)*m(0, 1)*m(1, 3) + m(2, 0)*m(0, 3)*m(1, 1);
        inv(3, 3) =  m(0, 0)*m(1, 1)*m(2, 2) - m(0, 0)*m(1, 2)*m(2, 1) - m(1, 0)*m(0, 1)*m(2, 2)
                     + m(1, 0)*m(0, 2)*m(2, 1) + m(2, 0)*m(0, 1)*m(1, 2) - m(2, 0)*m(0, 2)*m(1, 1);

        T det = m(0, 0)*inv(0, 0) + m(1, 0)*inv(1, 0) + m(2, 0)*inv(2, 0) + m(3, 0)*inv(3, 0);

        if (std::abs(det) < EPS) {
            std::cerr << "Matrix is not invertible!\n";
            *this = Matrix4x4<T>::create_identity();
            return;
        }

        T inv_det = 1.0 / det;
        for (int i = 0; i < 4; ++i)
            for (int j = 0; j < 4; ++j)
                (*this)(j, i) = inv(j, i) * inv_det;
    }

    // From ChatGPT lol
    void inverseTranspose() {
        Matrix4x4<T> temp = *this;

        // Step 1: Compute determinant via cofactor expansion along row 0
        T det = 0;
        for (int j = 0; j < 4; ++j) {
            Matrix<T, 3, 3> minor;
            for (int mi = 1; mi < 4; ++mi) {
                int r = mi - 1;
                int c = 0;
                for (int mj = 0; mj < 4; ++mj) {
                    if (mj == j) continue;
                    minor(c, r) = temp(mj, mi); // x=mj (col), y=mi (row)
                    c++;
                }
            }

            T cofactor =
                    minor(0, 0) * (minor(1, 1) * minor(2, 2) - minor(2, 1) * minor(1, 2)) -
                    minor(1, 0) * (minor(0, 1) * minor(2, 2) - minor(2, 1) * minor(0, 2)) +
                    minor(2, 0) * (minor(0, 1) * minor(1, 2) - minor(1, 1) * minor(0, 2));

            det += ((j % 2 == 0 ? 1 : -1) * temp(j, 0) * cofactor);
        }

        if (std::abs(det) < std::numeric_limits<T>::epsilon() * std::abs(det)) {
            std::cerr << "Matrix is not invertible!\n";
            *this = Matrix4x4<T>::create_identity();
            return;
        }

        // Step 2: Compute inverse and transpose at the same time (adjugate / det)
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                Matrix<T, 3, 3> minor;
                for (int mi = 0, r = 0; mi < 4; mi++) {
                    if (mi == i) continue;
                    for (int mj = 0, c = 0; mj < 4; mj++) {
                        if (mj == j) continue;
                        minor(c, r) = temp(mj, mi); // x=mj (col), y=mi (row)
                        c++;
                    }
                    r++;
                }

                T cofactor =
                        minor(0, 0) * (minor(1, 1) * minor(2, 2) - minor(2, 1) * minor(1, 2)) -
                        minor(1, 0) * (minor(0, 1) * minor(2, 2) - minor(2, 1) * minor(0, 2)) +
                        minor(2, 0) * (minor(0, 1) * minor(1, 2) - minor(1, 1) * minor(0, 2));

                (*this)(i, j) = ((i + j) % 2 == 0 ? 1 : -1) * cofactor / det; // Transposed: row i, col j
            }
        }
    }

    /* Operator overloads */
    /// Vector is always on the right
    Vec4<T> operator* (Vec4<T> const &other) {
        Vec4<T> result;
        for(int i = 0; i < 4; i++) {
            result[i] = other.x * this->data[i * 4 + 0] +
                        other.y * this->data[i * 4 + 1] +
                        other.z * this->data[i * 4 + 2] +
                        other.w * this->data[i * 4 + 3];
        }
        return result;
    }
    /// Direct to the multiply function
    inline Matrix4x4<T> operator*( Matrix4x4<T> const &other) { return this->multiply4x4(other); };

    template <typename U, size_t w, size_t h>
    friend Matrix<U, w, h> operator* ( Matrix4x4<U> const &cur,  Matrix<U, w, h> const &other);
};

template <typename T, size_t w, size_t h>
Matrix<T, w, h> operator* ( Matrix<T, w, h> const &cur,  Matrix4x4<T> const &other) {
 return cur.multiply(other);
}

template<typename T, size_t w, size_t h>
Matrix<T,w ,h> operator* ( Matrix4x4<T> const &cur,  Matrix<T, w, h> const &other) {
 return cur.multiply(other);
}

typedef Matrix4x4<float> Matrix4x4f;
static constexpr Matrix4x4f mat4x4_identity = {1., 0., 0., 0.,
                                           0., 1., 0., 0.,
                                           0., 0., 1., 0.,
                                           0., 0., 0., 1.};
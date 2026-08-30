//
// Created by iUV on 1/25/2026.
//
#include <algorithm>
#include <stdexcept>
#include "Math/geometry.h"
#include <iostream>

#pragma once
class RectBuffer {
public:
    RectBuffer() = delete;
    RectBuffer(int const w, const int h): width(w), size(w*h) {
        data = new float[size];
        resetBuffer();
    }
    ~RectBuffer() {
        delete[] data;
    }

    float &operator[] (size_t const i) {
//        if (i >= size) throw std::out_of_range("out of range!");
        //if (i >= size) std::cerr << "out of range!\n";
        size_t idx = std::min(i, size - 1);

        return data[idx];
    }

    void resetBuffer() {
        std::fill(data, data + size, -MAX_FLOAT); // set every value in zbuffer to -inf
    }
    size_t width;
    size_t size;
    float *data = nullptr;
};

//
// Created by iUV on 1/25/2026.
//
#include <algorithm>
#include <stdexcept>
#include "geometry.h"


#pragma once
class RectBuffer {
public:
    RectBuffer() { };
    RectBuffer(int const w, const int h): width(w), size(w*h) {
        data = new float[size];
        resetBuffer();
    }
    ~RectBuffer() {
        delete[] data;
    }

    float &operator[] (size_t const i) {
        if (i >= size) return data[size];
        return data[i];
    }

    void resetBuffer() {
        std::fill(data, data + size, -MAX_FLOAT); // set every value in zbuffer to -inf
    }
    size_t width;
    size_t size;
    float *data;
};

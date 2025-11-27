#pragma once

template <typename T>
struct Matrix {
    using value_type = T; 
    T* data;
    size_t rows, cols;
    
    // Constructor
    Matrix(size_t r, size_t c) : rows(r), cols(c) {
        data = new T[rows * cols];
    }
    
    // Destructor
    ~Matrix() {
        delete[] data;
    }

    // Move Constructor
    Matrix(Matrix&& other) noexcept 
        : data(other.data), rows(other.rows), cols(other.cols) {
        other.data = nullptr;
        other.rows = 0;
        other.cols = 0;
    }
    
    // Move Assignment Operator
    Matrix& operator=(Matrix&& other) noexcept {
        if (this != &other) {
            delete[] data;
            data = other.data;
            rows = other.rows;
            cols = other.cols;
            other.data = nullptr;
            other.rows = 0;
            other.cols = 0;
        }
        return *this;
    }
    
    // Access element (i,j)
    T& operator()(size_t i, size_t j) {
        return data[i * cols + j];
    }
    
    // Const access
    const T& operator()(size_t i, size_t j) const {
        return data[i * cols + j];
    }

    // Copy operations
    Matrix(const Matrix& other) : rows(other.rows), cols(other.cols) {
        data = new T[rows * cols];
        std::copy(other.data, other.data + rows * cols, data);
    }

    Matrix& operator=(const Matrix& other) {
        if (this != &other) {
            delete[] data;
            rows = other.rows;
            cols = other.cols;
            data = new T[rows * cols];
            std::copy(other.data, other.data + rows * cols, data);
        }
        return *this;
    }
};
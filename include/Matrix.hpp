#pragma once


template <typename T>
struct Matrix{
    T * data;
    int rows, cols;
    // Constructor
    Matrix(int r, int c) : rows(r), cols(c) {
        data = new T[rows * cols];
    }
    
    // Destructor
    ~Matrix() {
        delete[] data;
    }

    // Move Constructor
    Matrix(Matrix&& other) noexcept 
        : data(other.data), rows(other.rows), cols(other.cols) {
        other.data = nullptr;  // Important: leave source in valid state
        other.rows = 0;
        other.cols = 0;
    }
    
    // Move Assignment Operator
    Matrix& operator=(Matrix&& other) noexcept {
        if (this != &other) {
            // Clean up current resources
            delete[] data;
            
            // Transfer ownership
            data = other.data;
            rows = other.rows;
            cols = other.cols;
            
            // Leave source in valid state
            other.data = nullptr;
            other.rows = 0;
            other.cols = 0;
        }
        return *this;
    }
    
    // Access element (i,j)
    T& operator()(int i, int j) {
        return data[i * cols + j];
    }
    
    // Const access
    const T& operator()(int i, int j) const {
        return data[i * cols + j];
    }

    // Copy operations (if needed - but expensive!)
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


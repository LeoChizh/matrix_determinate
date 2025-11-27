#pragma once

template <typename T>
struct MatrixJagged {
    using value_type = T; 
    T** data;
    size_t rows, cols;
    
    // Constructor
    MatrixJagged(size_t r, size_t c) : rows(r), cols(c) {
        data = new T*[rows];
        for (size_t i = 0; i < rows; i++) {
            data[i] = new T[cols];
        }
    }
    
    // Destructor
    ~MatrixJagged() {
        clear();
    }
    
    // Copy Constructor
    MatrixJagged(const MatrixJagged& other) : rows(other.rows), cols(other.cols) {
        data = new T*[rows];
        for (size_t i = 0; i < rows; i++) {
            data[i] = new T[cols];
            std::copy(other.data[i], other.data[i] + cols, data[i]);
        }
    }
    
    // Copy Assignment Operator
    MatrixJagged& operator=(const MatrixJagged& other) {
        if (this != &other) {
            clear();
            rows = other.rows;
            cols = other.cols;
            data = new T*[rows];
            for (size_t i = 0; i < rows; i++) {
                data[i] = new T[cols];
                std::copy(other.data[i], other.data[i] + cols, data[i]);
            }
        }
        return *this;
    }
    
    // Move Constructor
    MatrixJagged(MatrixJagged&& other) noexcept 
        : data(other.data), rows(other.rows), cols(other.cols) {
        other.data = nullptr;
        other.rows = 0;
        other.cols = 0;
    }
    
    // Move Assignment Operator
    MatrixJagged& operator=(MatrixJagged&& other) noexcept {
        if (this != &other) {
            clear();
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
        return data[i][j];
    }
    
    // Const access
    const T& operator()(size_t i, size_t j) const {
        return data[i][j];
    }

private:
    void clear() {
        if (data) {
            for (size_t i = 0; i < rows; i++) {
                if (data[i]) {
                    delete[] data[i];
                }
            }
            delete[] data;
            data = nullptr;
        }
        rows = 0;
        cols = 0;
    }
};
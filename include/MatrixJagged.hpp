#pragma once

template <typename T>
struct MatrixJagged {
    T** data;
    int rows, cols;
    
    // Constructor
    MatrixJagged(int r, int c) : rows(r), cols(c) {
        data = new T*[rows];
        for (int i = 0; i < rows; i++) {
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
        for (int i = 0; i < rows; i++) {
            data[i] = new T[cols];
            std::copy(other.data[i], other.data[i] + cols, data[i]);
        }
    }
    
    // Copy Assignment Operator
    MatrixJagged& operator=(const MatrixJagged& other) {
        if (this != &other) {
            // Clean up current resources
            clear();
            
            // Copy dimensions
            rows = other.rows;
            cols = other.cols;
            
            // Allocate new memory and copy data
            data = new T*[rows];
            for (int i = 0; i < rows; i++) {
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
    T& operator()(int i, int j) {
        return data[i][j];
    }
    
    // Const access
    const T& operator()(int i, int j) const {
        return data[i][j];
    }

private:
    // Helper function to clean up resources
    void clear() {
        if (data) {
            for (int i = 0; i < rows; i++) {
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
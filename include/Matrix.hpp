#pragma once

template <typename T>
class Matrix {
private:
    T* data_;
    size_t rows_, cols_;

public:
    // Public read-only references to private members
    const size_t& rows;
    const size_t& cols;

    // Type alias needed for MatrixAlgorithms
    using value_type = T;

    // Constructor
    Matrix(size_t r, size_t c) : rows_(r), cols_(c), rows(rows_), cols(cols_) {
        size_t total_size = rows_ * cols_;
        if (total_size > 0) {
            data_ = new T[total_size];
            for (size_t i = 0; i < total_size; ++i) {
                data_[i] = T();
            }
        } else {
            data_ = nullptr;  // No allocation for zero-size matrix
        }
    }

    // Destructor
    ~Matrix() {
        delete[] data_;
    }

    // Move Constructor - must rebind references!
    Matrix(Matrix&& other) noexcept 
        : rows_(other.rows_), cols_(other.cols_), 
          rows(rows_), cols(cols_) {  // rebind to our own members
        data_ = other.data_;
        other.data_ = nullptr;
        other.rows_ = 0;
        other.cols_ = 0;
        // Note: other.rows and other.cols now reference other.rows_/cols_ which are 0
    }

    // Move Assignment - must rebind references!
    Matrix& operator=(Matrix&& other) noexcept {
        if (this != &other) {
            delete[] data_;
            
            // Update our private members
            rows_ = other.rows_;
            cols_ = other.cols_;
            data_ = other.data_;
            
            // References are already bound to our rows_/cols_ during construction
            // so they automatically reflect the new values
            
            other.data_ = nullptr;
            other.rows_ = 0;
            other.cols_ = 0;
        }
        return *this;
    }

    // Copy Constructor - must rebind references!
    Matrix(const Matrix& other) 
        : rows_(other.rows_), cols_(other.cols_), 
          rows(rows_), cols(cols_) {  // bind to our own members
        data_ = new T[rows_ * cols_];
        for (size_t i = 0; i < rows_ * cols_; ++i) {
            data_[i] = other.data_[i];
        }
    }

    // Copy Assignment - must rebind references!
    Matrix& operator=(const Matrix& other) {
        if (this != &other) {
            delete[] data_;
            
            // Update our private members  
            rows_ = other.rows_;
            cols_ = other.cols_;
            
            // References automatically reflect the new values
            
            data_ = new T[rows_ * cols_];
            for (size_t i = 0; i < rows_ * cols_; ++i) {
                data_[i] = other.data_[i];
            }
        }
        return *this;
    }

    // Accessors
    T& operator()(size_t row, size_t col) {
        return data_[row * cols_ + col];
    }
    
    const T& operator()(size_t row, size_t col) const {
        return data_[row * cols_ + col];
    }

    // Size getter
    size_t size() const { return rows_ * cols_; }

    // Additional useful methods
    void fill(const T& value) {
        for (size_t i = 0; i < rows_ * cols_; ++i) {
            data_[i] = value;
        }
    }

};
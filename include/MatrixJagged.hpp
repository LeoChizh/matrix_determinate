#pragma once

template <typename T>
class MatrixJagged {
private:
    T** data_;
    size_t rows_, cols_;

public:
    // Public read-only references to private members
    const size_t& rows;
    const size_t& cols;

    // Type alias needed for MatrixAlgorithms
    using value_type = T;

    // Constructor
    MatrixJagged(size_t r, size_t c) : rows_(r), cols_(c), rows(rows_), cols(cols_) {
        data_ = new T*[rows_];
        for (size_t i = 0; i < rows_; i++) {
            data_[i] = new T[cols_];
            for (size_t j = 0; j < cols_; j++) {
                data_[i][j] = T();
            }
        }
    }

    // Destructor
    ~MatrixJagged() {
        clear();
    }

    // Move Constructor - must rebind references!
    MatrixJagged(MatrixJagged&& other) noexcept 
        : rows_(other.rows_), cols_(other.cols_), 
          rows(rows_), cols(cols_) {  // rebind to our own members
        data_ = other.data_;
        other.data_ = nullptr;
        other.rows_ = 0;
        other.cols_ = 0;
        // Note: other.rows and other.cols now reference other.rows_/cols_ which are 0
    }

    // Move Assignment - must rebind references!
    MatrixJagged& operator=(MatrixJagged&& other) noexcept {
        if (this != &other) {
            clear();
            
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
    MatrixJagged(const MatrixJagged& other) 
        : rows_(other.rows_), cols_(other.cols_), 
          rows(rows_), cols(cols_) {  // bind to our own members
        data_ = new T*[rows_];
        for (size_t i = 0; i < rows_; i++) {
            data_[i] = new T[cols_];
            for (size_t j = 0; j < cols_; j++) {
                data_[i][j] = other.data_[i][j];
            }
        }
    }

    // Copy Assignment - must rebind references!
    MatrixJagged& operator=(const MatrixJagged& other) {
        if (this != &other) {
            clear();
            
            // Update our private members  
            rows_ = other.rows_;
            cols_ = other.cols_;
            
            // References automatically reflect the new values
            
            data_ = new T*[rows_];
            for (size_t i = 0; i < rows_; i++) {
                data_[i] = new T[cols_];
                for (size_t j = 0; j < cols_; j++) {
                    data_[i][j] = other.data_[i][j];
                }
            }
        }
        return *this;
    }

    // Accessors
    T& operator()(size_t row, size_t col) {
        return data_[row][col];
    }
    
    const T& operator()(size_t row, size_t col) const {
        return data_[row][col];
    }

    // Size getter
    size_t size() const { return rows_ * cols_; }

    // Additional useful methods
    void fill(const T& value) {
        for (size_t i = 0; i < rows_; ++i) {
            for (size_t j = 0; j < cols_; ++j) {
                data_[i][j] = value;
            }
        }
    }


private:
    void clear() {
        if (data_) {
            for (size_t i = 0; i < rows_; i++) {
                delete[] data_[i];
            }
            delete[] data_;
            data_ = nullptr;
        }
        rows_ = 0;
        cols_ = 0;
    }
};
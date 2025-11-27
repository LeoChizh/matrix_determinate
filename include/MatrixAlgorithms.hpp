#pragma once
#include <optional>
#include <vector>
#include <cmath>
#include <stdexcept>
#include <utility>
#include <cstdint>
#include <limits>
#include "Matrix.hpp"
#include "MatrixJagged.hpp"

template<typename MatrixType>
class MatrixAlgorithms {
public:
    using value_type = typename MatrixType::value_type;
    
    // Error codes for determinant calculation
    enum class DeterminantError {
        SUCCESS = 0,          // Successfully computed determinant
        NOT_SQUARE,           // Matrix is not square - cannot compute determinant
        OVERFLOW_NUMBER              // Determinant value overflows the return type
        // Note: Singular matrices are SUCCESS with value 0
    };
    
    // Result type for determinant calculation
    struct DeterminantResult {
        DeterminantError error;
        value_type value;
        
        bool success() const { return error == DeterminantError::SUCCESS; }
        operator bool() const { return success(); }
    };
    
    // Creates identity matrix of the same size as reference matrix
    // Returns std::nullopt if matrix is not square
    static std::optional<MatrixType> createIdentity(const MatrixType& reference) {
        if (reference.rows != reference.cols) {
            return std::nullopt;
        }
        
        MatrixType identity(reference.rows, reference.cols);
        
        for (size_t i = 0; i < identity.rows; ++i) {
            for (size_t j = 0; j < identity.cols; ++j) {
                identity(i, j) = (i == j) ? value_type{1} : value_type{0};
            }
        }
        
        return identity;
    }
    
private:
    // Convert any matrix to double matrix
    static Matrix<double> convertToDoubleMatrix(const MatrixType& matrix) {
        Matrix<double> result(matrix.rows, matrix.cols);
        for (size_t i = 0; i < matrix.rows; ++i) {
            for (size_t j = 0; j < matrix.cols; ++j) {
                result(i, j) = static_cast<double>(matrix(i, j));
            }
        }
        return result;
    }
    
    // Find the pivot element in a column for LU decomposition
    static size_t findPivot(const Matrix<double>& LU, size_t col, size_t start_row) {
        size_t max_index = start_row;
        double max_val = std::abs(LU(start_row, col));
        
        for (size_t i = start_row + 1; i < LU.rows; ++i) {
            double current_val = std::abs(LU(i, col));
            if (current_val > max_val) {
                max_val = current_val;
                max_index = i;
            }
        }
        
        return max_index;
    }
    
    // Swap two rows in the matrix and update pivot tracking
    static void swapRows(Matrix<double>& LU, size_t row1, size_t row2, 
                        std::vector<size_t>& pivot, int& sign) {
        if (row1 == row2) return;
        
        // Swap rows in LU matrix
        for (size_t j = 0; j < LU.cols; ++j) {
            std::swap(LU(row1, j), LU(row2, j));
        }
        
        // Update pivot vector and sign
        std::swap(pivot[row1], pivot[row2]);
        sign = -sign;
    }
    
    // Perform the elimination step for LU decomposition
    static void eliminateColumn(Matrix<double>& LU, size_t pivot_row) {
        const size_t n = LU.rows;
        double pivot_value = LU(pivot_row, pivot_row);
        
        for (size_t i = pivot_row + 1; i < n; ++i) {
            // Compute multiplier
            LU(i, pivot_row) /= pivot_value;
            
            // Update remaining submatrix
            for (size_t j = pivot_row + 1; j < n; ++j) {
                LU(i, j) -= LU(i, pivot_row) * LU(pivot_row, j);
            }
        }
    }
    
    // Calculate determinant from LU decomposition result with overflow check
    static DeterminantResult calculateDeterminantFromLUProtected(
        const Matrix<double>& LU, int sign) {
        
        double det_double = static_cast<double>(sign);
        const size_t n = LU.rows;
        
        for (size_t i = 0; i < n; ++i) {
            det_double *= LU(i, i);
        }
        
        // For integer types, check for overflow when converting back
        if constexpr (std::is_integral_v<value_type>) {
            constexpr double max_val = static_cast<double>(std::numeric_limits<value_type>::max());
            constexpr double min_val = static_cast<double>(std::numeric_limits<value_type>::min());
            
            if (det_double > max_val || det_double < min_val) {
                // Overflow detected
                return {DeterminantError::OVERFLOW_NUMBER, value_type{0}};
            }
            
            // Round to nearest integer for integral types
            return {DeterminantError::SUCCESS, static_cast<value_type>(std::round(det_double))};
        } else {
            // For floating-point types, direct conversion
            return {DeterminantError::SUCCESS, static_cast<value_type>(det_double)};
        }
    }
    
    // Initialize pivot vector
    static std::vector<size_t> initializePivot(size_t n) {
        std::vector<size_t> pivot(n);
        for (size_t i = 0; i < n; ++i) {
            pivot[i] = i;
        }
        return pivot;
    }
    
    // Check if matrix is singular - all potential pivots in column are zero
    static bool isColumnSingular(const Matrix<double>& LU, size_t col, size_t start_row) {
        constexpr double tolerance = 1e-12;
        for (size_t i = start_row; i < LU.rows; ++i) {
            if (std::abs(LU(i, col)) > tolerance) {
                return false; // Found a non-zero pivot candidate
            }
        }
        return true; // All zeros in this column
    }
    
public:
    // Calculate determinant using LU decomposition with partial pivoting
    static DeterminantResult determinantLU(const MatrixType& matrix) {
        // Validate input matrix
        if (matrix.rows != matrix.cols) {
            return {DeterminantError::NOT_SQUARE, value_type{0}};
        }
        
        const size_t n = matrix.rows;
        
        // Handle edge cases
        if (n == 0) return {DeterminantError::SUCCESS, value_type{1}};
        if (n == 1) return {DeterminantError::SUCCESS, matrix(0, 0)};
        
        // Normal LU decomposition (always uses double internally)
        Matrix<double> LU = convertToDoubleMatrix(matrix);
        std::vector<size_t> pivot = initializePivot(n);
        int sign = 1;
        
        // Perform LU decomposition with partial pivoting
        for (size_t k = 0; k < n; ++k) {
            // Check if the entire column is singular (all potential pivots are zero)
            if (isColumnSingular(LU, k, k)) {
                // Matrix is singular - determinant is 0 (SUCCESS case)
                return {DeterminantError::SUCCESS, value_type{0}};
            }
            
            // Find pivot element in current column
            size_t pivot_index = findPivot(LU, k, k);
            
            // Swap rows if necessary
            if (pivot_index != k) {
                swapRows(LU, k, pivot_index, pivot, sign);
            }
            
            // Perform elimination for current column
            eliminateColumn(LU, k);
        }
        
        // Use protected calculation with overflow check
        return calculateDeterminantFromLUProtected(LU, sign);
    }
};
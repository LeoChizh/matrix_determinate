#pragma once
#include <optional>
#include "Matrix.hpp"
#include "MatrixJagged.hpp"


template<typename MatrixType>
class MatrixAlgorithms {
public:
    using value_type = typename MatrixType::value_type;
    
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
};
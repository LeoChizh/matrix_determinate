#include <gtest/gtest.h>
#include <optional>
#include "MatrixAlgorithms.hpp"
#include "Matrix.hpp"
#include "MatrixJagged.hpp"

// Test fixture for MatrixAlgorithms
template<typename T>
class MatrixAlgorithmsTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Common setup if needed
    }

    void TearDown() override {
        // Common cleanup if needed
    }
};

// Test both Matrix and MatrixJagged types
using MatrixTypes = ::testing::Types<Matrix<int>, Matrix<double>, MatrixJagged<int>, MatrixJagged<float>>;
TYPED_TEST_SUITE(MatrixAlgorithmsTest, MatrixTypes);

// Test successful identity matrix creation for square matrices
TYPED_TEST(MatrixAlgorithmsTest, CreateIdentity_SquareMatrix_Success) {
    using MatrixType = TypeParam;
    using ValueType = typename MatrixType::value_type;
    
    // Test with different square sizes
    MatrixType square3x3(3, 3);
    MatrixType square1x1(1, 1);
    MatrixType square5x5(5, 5);
    
    // Create identity matrices
    auto identity3 = MatrixAlgorithms<MatrixType>::createIdentity(square3x3);
    auto identity1 = MatrixAlgorithms<MatrixType>::createIdentity(square1x1);
    auto identity5 = MatrixAlgorithms<MatrixType>::createIdentity(square5x5);
    
    // Verify results are present
    ASSERT_TRUE(identity3.has_value());
    ASSERT_TRUE(identity1.has_value());
    ASSERT_TRUE(identity5.has_value());
    
    // Verify dimensions
    EXPECT_EQ(identity3->rows, 3);
    EXPECT_EQ(identity3->cols, 3);
    EXPECT_EQ(identity1->rows, 1);
    EXPECT_EQ(identity1->cols, 1);
    EXPECT_EQ(identity5->rows, 5);
    EXPECT_EQ(identity5->cols, 5);
    
    // Verify identity matrix properties for 3x3
    for (size_t i = 0; i < 3; ++i) {
        for (size_t j = 0; j < 3; ++j) {
            ValueType expected = (i == j) ? ValueType{1} : ValueType{0};
            EXPECT_EQ((*identity3)(i, j), expected) << "Mismatch at (" << i << ", " << j << ")";
        }
    }
    
    // Verify 1x1 identity
    EXPECT_EQ((*identity1)(0, 0), ValueType{1});
    
    // Verify 5x5 identity diagonal
    for (size_t i = 0; i < 5; ++i) {
        EXPECT_EQ((*identity5)(i, i), ValueType{1}) << "Diagonal mismatch at (" << i << ", " << i << ")";
        for (size_t j = 0; j < 5; ++j) {
            if (i != j) {
                EXPECT_EQ((*identity5)(i, j), ValueType{0}) << "Off-diagonal not zero at (" << i << ", " << j << ")";
            }
        }
    }
}

// Test failure for non-square matrices
TYPED_TEST(MatrixAlgorithmsTest, CreateIdentity_NonSquareMatrix_ReturnsNullopt) {
    using MatrixType = TypeParam;
    
    // Test various non-square dimensions
    MatrixType rect2x3(2, 3);
    MatrixType rect3x2(3, 2);
    MatrixType rect1x5(1, 5);
    MatrixType rect4x1(4, 1);
    
    // All should return std::nullopt
    auto result1 = MatrixAlgorithms<MatrixType>::createIdentity(rect2x3);
    auto result2 = MatrixAlgorithms<MatrixType>::createIdentity(rect3x2);
    auto result3 = MatrixAlgorithms<MatrixType>::createIdentity(rect1x5);
    auto result4 = MatrixAlgorithms<MatrixType>::createIdentity(rect4x1);
    
    EXPECT_FALSE(result1.has_value());
    EXPECT_FALSE(result2.has_value());
    EXPECT_FALSE(result3.has_value());
    EXPECT_FALSE(result4.has_value());
}

// Test that original matrix data is not modified
TYPED_TEST(MatrixAlgorithmsTest, CreateIdentity_OriginalMatrixNotModified) {
    using MatrixType = TypeParam;
    using ValueType = typename MatrixType::value_type;
    
    // Create a square matrix and fill it with specific values
    MatrixType original(2, 2);
    original(0, 0) = ValueType{10};
    original(0, 1) = ValueType{20};
    original(1, 0) = ValueType{30};
    original(1, 1) = ValueType{40};
    
    // Store original values for comparison
    ValueType original_00 = original(0, 0);
    ValueType original_01 = original(0, 1);
    ValueType original_10 = original(1, 0);
    ValueType original_11 = original(1, 1);
    
    // Create identity matrix
    auto identity = MatrixAlgorithms<MatrixType>::createIdentity(original);
    
    // Verify original matrix is unchanged
    EXPECT_EQ(original(0, 0), original_00);
    EXPECT_EQ(original(0, 1), original_01);
    EXPECT_EQ(original(1, 0), original_10);
    EXPECT_EQ(original(1, 1), original_11);
    
    // Verify identity matrix is correct (and different from original)
    ASSERT_TRUE(identity.has_value());
    EXPECT_EQ((*identity)(0, 0), ValueType{1});
    EXPECT_EQ((*identity)(0, 1), ValueType{0});
    EXPECT_EQ((*identity)(1, 0), ValueType{0});
    EXPECT_EQ((*identity)(1, 1), ValueType{1});
}

// Test with zero-sized matrix (edge case)
TYPED_TEST(MatrixAlgorithmsTest, CreateIdentity_ZeroSizedMatrix_Success) {
    using MatrixType = TypeParam;
    
    MatrixType zeroMatrix(0, 0);
    
    auto identity = MatrixAlgorithms<MatrixType>::createIdentity(zeroMatrix);
    
    ASSERT_TRUE(identity.has_value());
    EXPECT_EQ(identity->rows, 0);
    EXPECT_EQ(identity->cols, 0);
}

// Test that returned matrix is independent (deep copy)
TYPED_TEST(MatrixAlgorithmsTest, CreateIdentity_ReturnedMatrixIsIndependent) {
    using MatrixType = TypeParam;
    using ValueType = typename MatrixType::value_type;
    
    MatrixType original(2, 2);
    original(0, 0) = ValueType{5};
    original(1, 1) = ValueType{5};
    
    auto identity = MatrixAlgorithms<MatrixType>::createIdentity(original);
    ASSERT_TRUE(identity.has_value());
    
    // Modify the identity matrix
    (*identity)(0, 0) = ValueType{99};
    (*identity)(1, 0) = ValueType{99};
    
    // Verify original is unchanged
    EXPECT_EQ(original(0, 0), ValueType{5});
    EXPECT_EQ(original(1, 1), ValueType{5});
}

// Specific tests for different numeric types
TEST(MatrixAlgorithmsSpecific, DifferentNumericTypes) {
    // Test with integer matrix
    Matrix<int> intMatrix(2, 2);
    auto intIdentity = MatrixAlgorithms<Matrix<int>>::createIdentity(intMatrix);
    ASSERT_TRUE(intIdentity.has_value());
    EXPECT_EQ((*intIdentity)(0, 0), 1);
    EXPECT_EQ((*intIdentity)(1, 1), 1);
    
    // Test with double matrix
    Matrix<double> doubleMatrix(2, 2);
    auto doubleIdentity = MatrixAlgorithms<Matrix<double>>::createIdentity(doubleMatrix);
    ASSERT_TRUE(doubleIdentity.has_value());
    EXPECT_DOUBLE_EQ((*doubleIdentity)(0, 0), 1.0);
    EXPECT_DOUBLE_EQ((*doubleIdentity)(1, 1), 1.0);
    
    // Test with float jagged matrix
    MatrixJagged<float> floatJagged(2, 2);
    auto floatIdentity = MatrixAlgorithms<MatrixJagged<float>>::createIdentity(floatJagged);
    ASSERT_TRUE(floatIdentity.has_value());
    EXPECT_FLOAT_EQ((*floatIdentity)(0, 0), 1.0f);
    EXPECT_FLOAT_EQ((*floatIdentity)(1, 1), 1.0f);
}

// Test performance with large matrix
TYPED_TEST(MatrixAlgorithmsTest, CreateIdentity_LargeMatrix_Success) {
    using MatrixType = TypeParam;
    
    const size_t largeSize = 100;
    MatrixType largeMatrix(largeSize, largeSize);
    
    auto identity = MatrixAlgorithms<MatrixType>::createIdentity(largeMatrix);
    
    ASSERT_TRUE(identity.has_value());
    EXPECT_EQ(identity->rows, largeSize);
    EXPECT_EQ(identity->cols, largeSize);
    
    // Verify a few key elements
    EXPECT_EQ((*identity)(0, 0), typename MatrixType::value_type{1});
    EXPECT_EQ((*identity)(largeSize-1, largeSize-1), typename MatrixType::value_type{1});
    EXPECT_EQ((*identity)(0, 1), typename MatrixType::value_type{0});
    EXPECT_EQ((*identity)(1, 0), typename MatrixType::value_type{0});
}

// Test determinant calculation using LU decomposition
TYPED_TEST(MatrixAlgorithmsTest, DeterminantLU_2x2IntegerMatrix_Success) {
    using MatrixType = TypeParam;
    using ValueType = typename MatrixType::value_type;
    
    MatrixType mat(2, 2);
    mat(0, 0) = ValueType{1}; mat(0, 1) = ValueType{2};
    mat(1, 0) = ValueType{3}; mat(1, 1) = ValueType{4};
    
    // det = 1*4 - 2*3 = -2
    auto result = MatrixAlgorithms<MatrixType>::determinantLU(mat);
    
    EXPECT_TRUE(result.first); // Success
    EXPECT_NEAR(result.second, -2.0, 1e-10);
}

TYPED_TEST(MatrixAlgorithmsTest, DeterminantLU_NonSquareMatrix_Failure) {
    using MatrixType = TypeParam;
    
    MatrixType rect2x3(2, 3);
    
    auto result = MatrixAlgorithms<MatrixType>::determinantLU(rect2x3);
    
    EXPECT_FALSE(result.first); // Failure
    EXPECT_EQ(result.second, 0.0); // Default value on failure
}

TYPED_TEST(MatrixAlgorithmsTest, DeterminantLU_SingularMatrix_Failure) {
    using MatrixType = TypeParam;
    using ValueType = typename MatrixType::value_type;
    
    MatrixType mat(3, 3);
    // Create a singular matrix (rows 1 and 2 are linearly dependent)
    mat(0, 0) = ValueType{1}; mat(0, 1) = ValueType{2}; mat(0, 2) = ValueType{3};
    mat(1, 0) = ValueType{4}; mat(1, 1) = ValueType{5}; mat(1, 2) = ValueType{6};
    mat(2, 0) = ValueType{7}; mat(2, 1) = ValueType{8}; mat(2, 2) = ValueType{9};
    
    auto result = MatrixAlgorithms<MatrixType>::determinantLU(mat);
    
    EXPECT_TRUE(result.first); // Failure - matrix is singular
    EXPECT_EQ(result.second, 0.0); // Default value on failure
}

TYPED_TEST(MatrixAlgorithmsTest, DeterminantLU_3x3IntegerMatrix_Success) {
    using MatrixType = TypeParam;
    using ValueType = typename MatrixType::value_type;
    
    MatrixType mat(3, 3);
    mat(0, 0) = ValueType{2}; mat(0, 1) = ValueType{-1}; mat(0, 2) = ValueType{0};
    mat(1, 0) = ValueType{-1}; mat(1, 1) = ValueType{2}; mat(1, 2) = ValueType{-1};
    mat(2, 0) = ValueType{0}; mat(2, 1) = ValueType{-1}; mat(2, 2) = ValueType{2};
    
    // det = 2*(2*2 - (-1)*(-1)) - (-1)*(-1*2 - (-1)*0) + 0*... = 2*(4-1) - (-1)*(-2) = 6 - 2 = 4
    auto result = MatrixAlgorithms<MatrixType>::determinantLU(mat);
    
    EXPECT_TRUE(result.first);
    EXPECT_NEAR(result.second, 4.0, 1e-10);
}

TYPED_TEST(MatrixAlgorithmsTest, DeterminantLU_IdentityMatrix_ReturnsOne) {
    using MatrixType = TypeParam;
    
    MatrixType identity(3, 3);
    auto identity_opt = MatrixAlgorithms<MatrixType>::createIdentity(identity);
    ASSERT_TRUE(identity_opt.has_value());
    
    auto result = MatrixAlgorithms<MatrixType>::determinantLU(*identity_opt);
    
    EXPECT_TRUE(result.first);
    EXPECT_NEAR(result.second, 1.0, 1e-10);
}

TYPED_TEST(MatrixAlgorithmsTest, DeterminantLU_ZeroSizedMatrix_ReturnsOne) {
    using MatrixType = TypeParam;
    
    MatrixType zeroMatrix(0, 0);
    
    auto result = MatrixAlgorithms<MatrixType>::determinantLU(zeroMatrix);
    
    EXPECT_TRUE(result.first);
    EXPECT_NEAR(result.second, 1.0, 1e-10);
}

TYPED_TEST(MatrixAlgorithmsTest, DeterminantLU_1x1IntegerMatrix_Success) {
    using MatrixType = TypeParam;
    using ValueType = typename MatrixType::value_type;
    
    MatrixType mat(1, 1);
    mat(0, 0) = ValueType{5};
    
    auto result = MatrixAlgorithms<MatrixType>::determinantLU(mat);
    
    EXPECT_TRUE(result.first);
    EXPECT_NEAR(result.second, 5.0, 1e-10);
}

TYPED_TEST(MatrixAlgorithmsTest, DeterminantLU_WithRowSwaps_Success) {
    using MatrixType = TypeParam;
    using ValueType = typename MatrixType::value_type;
    
    MatrixType mat(3, 3);
    // Matrix that requires row swapping during LU decomposition
    mat(0, 0) = ValueType{0}; mat(0, 1) = ValueType{1}; mat(0, 2) = ValueType{0};
    mat(1, 0) = ValueType{1}; mat(1, 1) = ValueType{0}; mat(1, 2) = ValueType{0};
    mat(2, 0) = ValueType{0}; mat(2, 1) = ValueType{0}; mat(2, 2) = ValueType{1};
    
    // det = -1 (due to row swap)
    auto result = MatrixAlgorithms<MatrixType>::determinantLU(mat);
    
    EXPECT_TRUE(result.first);
    EXPECT_NEAR(result.second, -1.0, 1e-10);
}

TYPED_TEST(MatrixAlgorithmsTest, DeterminantLU_LargeIntegerMatrix_Success) {
    using MatrixType = TypeParam;
    using ValueType = typename MatrixType::value_type;
    
    const size_t size = 5;
    MatrixType mat(size, size);
    
    // Create a diagonally dominant matrix (non-singular)
    for (size_t i = 0; i < size; ++i) {
        for (size_t j = 0; j < size; ++j) {
            if (i == j) {
                mat(i, j) = ValueType{size + 1};
            } else {
                mat(i, j) = ValueType{1};
            }
        }
    }
    
    auto result = MatrixAlgorithms<MatrixType>::determinantLU(mat);
    
    EXPECT_TRUE(result.first);
    // For diagonally dominant matrix, determinant should be non-zero
    EXPECT_GT(std::abs(result.second), 0.0);
}

// Test that original matrix is not modified
TYPED_TEST(MatrixAlgorithmsTest, DeterminantLU_OriginalMatrixNotModified) {
    using MatrixType = TypeParam;
    using ValueType = typename MatrixType::value_type;
    
    MatrixType original(2, 2);
    original(0, 0) = ValueType{1}; original(0, 1) = ValueType{2};
    original(1, 0) = ValueType{3}; original(1, 1) = ValueType{4};
    
    // Store original values
    ValueType original_00 = original(0, 0);
    ValueType original_01 = original(0, 1);
    ValueType original_10 = original(1, 0);
    ValueType original_11 = original(1, 1);
    
    auto result = MatrixAlgorithms<MatrixType>::determinantLU(original);
    
    // Verify original matrix is unchanged
    EXPECT_EQ(original(0, 0), original_00);
    EXPECT_EQ(original(0, 1), original_01);
    EXPECT_EQ(original(1, 0), original_10);
    EXPECT_EQ(original(1, 1), original_11);
    
    // Verify determinant is correct
    EXPECT_TRUE(result.first);
    EXPECT_NEAR(result.second, -2.0, 1e-10);
}
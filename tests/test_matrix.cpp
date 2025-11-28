#include <gtest/gtest.h>
#include "Matrix.hpp"  // Your Matrix header file
#include "MatrixJagged.hpp"  // Your MatrixJagged header file

// Test fixture template to test both Matrix and MatrixJagged
template<typename T>
class MatrixTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Common setup for tests
    }

    void TearDown() override {
        // Common cleanup for tests
    }
};

// Define the types we want to test
using MatrixTypes = ::testing::Types<Matrix<int>, MatrixJagged<int>>;
TYPED_TEST_SUITE(MatrixTest, MatrixTypes);

// Test basic construction and dimensions
TYPED_TEST(MatrixTest, ConstructorAndDimensions) {
    TypeParam mat(3, 4);
    
    EXPECT_EQ(mat.rows, 3);
    EXPECT_EQ(mat.cols, 4);
    EXPECT_EQ(mat.size(), 12);  // Use size() instead of accessing data directly
}

// Test element access and modification
TYPED_TEST(MatrixTest, ElementAccess) {
    TypeParam mat(2, 2);
    
    // Test write access
    mat(0, 0) = 1;
    mat(0, 1) = 2;
    mat(1, 0) = 3;
    mat(1, 1) = 4;
    
    // Test read access
    EXPECT_EQ(mat(0, 0), 1);
    EXPECT_EQ(mat(0, 1), 2);
    EXPECT_EQ(mat(1, 0), 3);
    EXPECT_EQ(mat(1, 1), 4);
    
    // Test const access
    const TypeParam& const_mat = mat;
    EXPECT_EQ(const_mat(0, 0), 1);
}

// Test copy constructor
TYPED_TEST(MatrixTest, CopyConstructor) {
    TypeParam original(2, 2);
    original(0, 0) = 5;
    original(0, 1) = 6;
    original(1, 0) = 7;
    original(1, 1) = 8;
    
    // Create copy
    TypeParam copy = original;
    
    // Verify dimensions and data are copied
    EXPECT_EQ(copy.rows, 2);
    EXPECT_EQ(copy.cols, 2);
    EXPECT_EQ(copy(0, 0), 5);
    EXPECT_EQ(copy(0, 1), 6);
    EXPECT_EQ(copy(1, 0), 7);
    EXPECT_EQ(copy(1, 1), 8);
    
    // Verify deep copy - modifying copy doesn't affect original
    copy(0, 0) = 100;
    EXPECT_EQ(original(0, 0), 5);  // Original unchanged
    EXPECT_EQ(copy(0, 0), 100);    // Copy changed
}

// Test copy assignment operator
TYPED_TEST(MatrixTest, CopyAssignment) {
    TypeParam original(2, 2);
    original(0, 0) = 10;
    original(0, 1) = 20;
    
    TypeParam assigned(1, 1);  // Different size initially
    assigned(0, 0) = 999;
    
    // Perform copy assignment
    assigned = original;
    
    // Verify dimensions and data are copied
    EXPECT_EQ(assigned.rows, 2);
    EXPECT_EQ(assigned.cols, 2);
    EXPECT_EQ(assigned(0, 0), 10);
    EXPECT_EQ(assigned(0, 1), 20);
    
    // Verify deep copy
    assigned(0, 0) = 200;
    EXPECT_EQ(original(0, 0), 10);  // Original unchanged
    
    // Test self-assignment
    assigned = assigned;
    EXPECT_EQ(assigned(0, 0), 200);  // Should remain unchanged
}

// Test move constructor
TYPED_TEST(MatrixTest, MoveConstructor) {
    TypeParam original(2, 2);
    original(0, 0) = 30;
    original(0, 1) = 40;
    size_t original_rows = original.rows;
    size_t original_cols = original.cols;
    
    // Move construct
    TypeParam moved = std::move(original);
    
    // Verify moved object has the data
    EXPECT_EQ(moved.rows, original_rows);
    EXPECT_EQ(moved.cols, original_cols);
    EXPECT_EQ(moved(0, 0), 30);
    EXPECT_EQ(moved(0, 1), 40);
    
    // Verify original is in valid but empty state
    EXPECT_EQ(original.rows, 0);
    EXPECT_EQ(original.cols, 0);
}

// Test move assignment operator
TYPED_TEST(MatrixTest, MoveAssignment) {
    TypeParam original(2, 2);
    original(0, 0) = 50;
    original(0, 1) = 60;
    
    TypeParam assigned(1, 1);
    assigned(0, 0) = 999;
    
    // Move assign
    assigned = std::move(original);
    
    // Verify assigned object has the data
    EXPECT_EQ(assigned.rows, 2);
    EXPECT_EQ(assigned.cols, 2);
    EXPECT_EQ(assigned(0, 0), 50);
    EXPECT_EQ(assigned(0, 1), 60);
    
    // Verify original is in valid but empty state
    EXPECT_EQ(original.rows, 0);
    EXPECT_EQ(original.cols, 0);
}

// Test destructor doesn't cause double deletion
TYPED_TEST(MatrixTest, DestructorSafety) {
    // This test verifies that moved-from objects can be safely destroyed
    TypeParam* original = new TypeParam(2, 2);
    (*original)(0, 0) = 100;
    
    TypeParam moved = std::move(*original);
    delete original;  // Should not cause double deletion
    
    EXPECT_EQ(moved(0, 0), 100);
}

// Test with different data types
TYPED_TEST(MatrixTest, DifferentDataTypes) {
    // Note: We need to create a new type alias for different value types
    // This is a bit tricky with typed tests, so we'll test double separately
}

// Test large matrix operations
TYPED_TEST(MatrixTest, LargeMatrix) {
    const size_t large_size = 100;
    TypeParam large_mat(large_size, large_size);
    
    // Initialize large matrix
    for (size_t i = 0; i < large_size; i++) {
        for (size_t j = 0; j < large_size; j++) {
            large_mat(i, j) = static_cast<int>(i * large_size + j);
        }
    }
    
    // Verify some values
    EXPECT_EQ(large_mat(0, 0), 0);
    EXPECT_EQ(large_mat(large_size-1, large_size-1), static_cast<int>(large_size*large_size - 1));
    
    // Test copy of large matrix
    TypeParam large_copy = large_mat;
    EXPECT_EQ(large_copy(large_size-1, large_size-1), static_cast<int>(large_size*large_size - 1));
    
    // Test move of large matrix
    TypeParam large_moved = std::move(large_mat);
    EXPECT_EQ(large_moved(large_size-1, large_size-1), static_cast<int>(large_size*large_size - 1));
}

// Test edge cases
TYPED_TEST(MatrixTest, EdgeCases) {
    // 1x1 matrix
    TypeParam single(1, 1);
    single(0, 0) = 42;
    EXPECT_EQ(single(0, 0), 42);
    
    // Test fill method
    TypeParam filled(3, 3);
    filled.fill(7);
    for (size_t i = 0; i < 3; i++) {
        for (size_t j = 0; j < 3; j++) {
            EXPECT_EQ(filled(i, j), 7);
        }
    }
}

// Additional tests for specific Matrix functionality
TEST(MatrixSpecificTest, MemoryLayout) {
    // This test only makes sense for the contiguous Matrix class
    Matrix<int> mat(2, 3);
    
    // Fill matrix
    mat(0, 0) = 1; mat(0, 1) = 2; mat(0, 2) = 3;
    mat(1, 0) = 4; mat(1, 1) = 5; mat(1, 2) = 6;
    
    // For MatrixJagged, we can't test the internal memory layout directly
    // since data is private and the layout is different
}

// Test double type specifically
TEST(MatrixDoubleTest, DoubleOperations) {
    Matrix<double> double_mat(2, 2);
    double_mat(0, 0) = 1.5;
    double_mat(0, 1) = 2.5;
    
    EXPECT_DOUBLE_EQ(double_mat(0, 0), 1.5);
    EXPECT_DOUBLE_EQ(double_mat(0, 1), 2.5);
}

TEST(MatrixJaggedDoubleTest, DoubleOperations) {
    MatrixJagged<double> double_mat(2, 2);
    double_mat(0, 0) = 1.5;
    double_mat(0, 1) = 2.5;
    
    EXPECT_DOUBLE_EQ(double_mat(0, 0), 1.5);
    EXPECT_DOUBLE_EQ(double_mat(0, 1), 2.5);
}

// Test very large matrices (stress test)
TEST(MatrixStressTest, VeryLargeMatrix) {
    Matrix<int> large(1000, 1000);
    large(999, 999) = 123;
    EXPECT_EQ(large(999, 999), 123);
}

TEST(MatrixJaggedStressTest, VeryLargeMatrix) {
    MatrixJagged<int> large(1000, 1000);
    large(999, 999) = 123;
    EXPECT_EQ(large(999, 999), 123);
}
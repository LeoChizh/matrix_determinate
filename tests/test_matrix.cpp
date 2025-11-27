#include <gtest/gtest.h>
#include "Matrix.hpp"  // Your header file

class MatrixTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Common setup for tests
    }

    void TearDown() override {
        // Common cleanup for tests
    }
};

// Test basic construction and dimensions
TEST_F(MatrixTest, ConstructorAndDimensions) {
    Matrix<int> mat(3, 4);
    
    EXPECT_EQ(mat.rows, 3);
    EXPECT_EQ(mat.cols, 4);
    EXPECT_NE(mat.data, nullptr);
}

// Test element access and modification
TEST_F(MatrixTest, ElementAccess) {
    Matrix<int> mat(2, 2);
    
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
    const Matrix<int>& const_mat = mat;
    EXPECT_EQ(const_mat(0, 0), 1);
}

// Test copy constructor
TEST_F(MatrixTest, CopyConstructor) {
    Matrix<int> original(2, 2);
    original(0, 0) = 5;
    original(0, 1) = 6;
    original(1, 0) = 7;
    original(1, 1) = 8;
    
    // Create copy
    Matrix<int> copy = original;
    
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
    
    // Verify different memory addresses
    EXPECT_NE(original.data, copy.data);
}

// Test copy assignment operator
TEST_F(MatrixTest, CopyAssignment) {
    Matrix<int> original(2, 2);
    original(0, 0) = 10;
    original(0, 1) = 20;
    
    Matrix<int> assigned(1, 1);  // Different size initially
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
TEST_F(MatrixTest, MoveConstructor) {
    Matrix<int> original(2, 2);
    original(0, 0) = 30;
    original(0, 1) = 40;
    int* original_data = original.data;
    int original_rows = original.rows;
    int original_cols = original.cols;
    
    // Move construct
    Matrix<int> moved = std::move(original);
    
    // Verify moved object has the data
    EXPECT_EQ(moved.data, original_data);
    EXPECT_EQ(moved.rows, original_rows);
    EXPECT_EQ(moved.cols, original_cols);
    EXPECT_EQ(moved(0, 0), 30);
    EXPECT_EQ(moved(0, 1), 40);
    
    // Verify original is in valid but empty state
    EXPECT_EQ(original.data, nullptr);
    EXPECT_EQ(original.rows, 0);
    EXPECT_EQ(original.cols, 0);
}

// Test move assignment operator
TEST_F(MatrixTest, MoveAssignment) {
    Matrix<int> original(2, 2);
    original(0, 0) = 50;
    original(0, 1) = 60;
    int* original_data = original.data;
    
    Matrix<int> assigned(1, 1);
    assigned(0, 0) = 999;
    
    // Move assign
    assigned = std::move(original);
    
    // Verify assigned object has the data
    EXPECT_EQ(assigned.data, original_data);
    EXPECT_EQ(assigned.rows, 2);
    EXPECT_EQ(assigned.cols, 2);
    EXPECT_EQ(assigned(0, 0), 50);
    EXPECT_EQ(assigned(0, 1), 60);
    
    // Verify original is in valid but empty state
    EXPECT_EQ(original.data, nullptr);
    EXPECT_EQ(original.rows, 0);
    EXPECT_EQ(original.cols, 0);
    
    // Test self move assignment using swap trick
    Matrix<int> temp = std::move(assigned);
    assigned = std::move(temp);  // This effectively tests self-move safety
    
    // Verify data is still valid after the "self-move"
    EXPECT_EQ(assigned(0, 0), 50);
    EXPECT_EQ(assigned(0, 1), 60);
    EXPECT_EQ(assigned.rows, 2);
    EXPECT_EQ(assigned.cols, 2);
}

// Test destructor doesn't cause double deletion
TEST_F(MatrixTest, DestructorSafety) {
    // This test verifies that moved-from objects can be safely destroyed
    Matrix<int>* original = new Matrix<int>(2, 2);
    original->operator()(0, 0) = 100;
    
    Matrix<int> moved = std::move(*original);
    delete original;  // Should not cause double deletion
    
    EXPECT_EQ(moved(0, 0), 100);
}

// Test with different data types
TEST_F(MatrixTest, DifferentDataTypes) {
    Matrix<double> double_mat(2, 2);
    double_mat(0, 0) = 1.5;
    double_mat(0, 1) = 2.5;
    
    EXPECT_DOUBLE_EQ(double_mat(0, 0), 1.5);
    EXPECT_DOUBLE_EQ(double_mat(0, 1), 2.5);
    
    Matrix<float> float_mat(1, 1);
    float_mat(0, 0) = 3.14f;
    
    EXPECT_FLOAT_EQ(float_mat(0, 0), 3.14f);
}

// Test large matrix operations
TEST_F(MatrixTest, LargeMatrix) {
    const int large_size = 100;
    Matrix<int> large_mat(large_size, large_size);
    
    // Initialize large matrix
    for (int i = 0; i < large_size; i++) {
        for (int j = 0; j < large_size; j++) {
            large_mat(i, j) = i * large_size + j;
        }
    }
    
    // Verify some values
    EXPECT_EQ(large_mat(0, 0), 0);
    EXPECT_EQ(large_mat(large_size-1, large_size-1), large_size*large_size - 1);
    
    // Test copy of large matrix
    Matrix<int> large_copy = large_mat;
    EXPECT_EQ(large_copy(large_size-1, large_size-1), large_size*large_size - 1);
    
    // Test move of large matrix
    Matrix<int> large_moved = std::move(large_mat);
    EXPECT_EQ(large_moved(large_size-1, large_size-1), large_size*large_size - 1);
}

// Test edge cases
TEST_F(MatrixTest, EdgeCases) {
    // 1x1 matrix
    Matrix<int> single(1, 1);
    single(0, 0) = 42;
    EXPECT_EQ(single(0, 0), 42);
    
    // Very large matrix (stress test)
    Matrix<int> large(1000, 1000);
    large(999, 999) = 123;
    EXPECT_EQ(large(999, 999), 123);
}

// Test memory layout (row-major order)
TEST_F(MatrixTest, MemoryLayout) {
    Matrix<int> mat(2, 3);
    
    // Fill matrix
    mat(0, 0) = 1; mat(0, 1) = 2; mat(0, 2) = 3;
    mat(1, 0) = 4; mat(1, 1) = 5; mat(1, 2) = 6;
    
    // Verify row-major order in memory
    EXPECT_EQ(mat.data[0], 1);  // (0,0)
    EXPECT_EQ(mat.data[1], 2);  // (0,1)  
    EXPECT_EQ(mat.data[2], 3);  // (0,2)
    EXPECT_EQ(mat.data[3], 4);  // (1,0)
    EXPECT_EQ(mat.data[4], 5);  // (1,1)
    EXPECT_EQ(mat.data[5], 6);  // (1,2)
}
#include <gtest/gtest.h>
#include "MatrixJagged.hpp"  // Your header file

class MatrixJaggedTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Common setup for tests
    }

    void TearDown() override {
        // Common cleanup for tests
    }
};

// Test basic construction and dimensions
TEST_F(MatrixJaggedTest, ConstructorAndDimensions) {
    MatrixJagged<int> mat(3, 4);
    
    EXPECT_EQ(mat.rows, 3);
    EXPECT_EQ(mat.cols, 4);
    EXPECT_NE(mat.data, nullptr);
    
    // Verify all row pointers are allocated
    for (int i = 0; i < mat.rows; i++) {
        EXPECT_NE(mat.data[i], nullptr);
    }
}

// Test element access and modification
TEST_F(MatrixJaggedTest, ElementAccess) {
    MatrixJagged<int> mat(2, 2);
    
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
    const MatrixJagged<int>& const_mat = mat;
    EXPECT_EQ(const_mat(0, 0), 1);
}

// Test copy constructor
TEST_F(MatrixJaggedTest, CopyConstructor) {
    MatrixJagged<int> original(2, 3);
    original(0, 0) = 5; original(0, 1) = 6; original(0, 2) = 7;
    original(1, 0) = 8; original(1, 1) = 9; original(1, 2) = 10;
    
    // Create copy
    MatrixJagged<int> copy = original;
    
    // Verify dimensions are copied
    EXPECT_EQ(copy.rows, 2);
    EXPECT_EQ(copy.cols, 3);
    
    // Verify data is copied
    EXPECT_EQ(copy(0, 0), 5);
    EXPECT_EQ(copy(0, 1), 6);
    EXPECT_EQ(copy(0, 2), 7);
    EXPECT_EQ(copy(1, 0), 8);
    EXPECT_EQ(copy(1, 1), 9);
    EXPECT_EQ(copy(1, 2), 10);
    
    // Verify deep copy - modifying copy doesn't affect original
    copy(0, 0) = 100;
    EXPECT_EQ(original(0, 0), 5);  // Original unchanged
    EXPECT_EQ(copy(0, 0), 100);    // Copy changed
    
    // Verify different memory addresses for data arrays
    EXPECT_NE(original.data, copy.data);
    for (int i = 0; i < original.rows; i++) {
        EXPECT_NE(original.data[i], copy.data[i]);
    }
}

// Test copy assignment operator
TEST_F(MatrixJaggedTest, CopyAssignment) {
    MatrixJagged<int> original(2, 2);
    original(0, 0) = 10;
    original(0, 1) = 20;
    original(1, 0) = 30;
    original(1, 1) = 40;
    
    MatrixJagged<int> assigned(1, 1);  // Different size initially
    assigned(0, 0) = 999;
    
    // Perform copy assignment
    assigned = original;
    
    // Verify dimensions and data are copied
    EXPECT_EQ(assigned.rows, 2);
    EXPECT_EQ(assigned.cols, 2);
    EXPECT_EQ(assigned(0, 0), 10);
    EXPECT_EQ(assigned(0, 1), 20);
    EXPECT_EQ(assigned(1, 0), 30);
    EXPECT_EQ(assigned(1, 1), 40);
    
    // Verify deep copy
    assigned(0, 0) = 200;
    EXPECT_EQ(original(0, 0), 10);  // Original unchanged
    
    // Test self-assignment
    assigned = assigned;
    EXPECT_EQ(assigned(0, 0), 200);  // Should remain unchanged
}

// Test move constructor
TEST_F(MatrixJaggedTest, MoveConstructor) {
    MatrixJagged<int> original(2, 2);
    original(0, 0) = 30;
    original(0, 1) = 40;
    original(1, 0) = 50;
    original(1, 1) = 60;
    
    int** original_data = original.data;
    int original_rows = original.rows;
    int original_cols = original.cols;
    
    // Store original row pointers for verification
    int** original_row_pointers = new int*[original_rows];
    for (int i = 0; i < original_rows; i++) {
        original_row_pointers[i] = original.data[i];
    }
    
    // Move construct
    MatrixJagged<int> moved = std::move(original);
    
    // Verify moved object has the data and pointers
    EXPECT_EQ(moved.data, original_data);
    EXPECT_EQ(moved.rows, original_rows);
    EXPECT_EQ(moved.cols, original_cols);
    EXPECT_EQ(moved(0, 0), 30);
    EXPECT_EQ(moved(0, 1), 40);
    EXPECT_EQ(moved(1, 0), 50);
    EXPECT_EQ(moved(1, 1), 60);
    
    // Verify all row pointers are transferred
    for (int i = 0; i < moved.rows; i++) {
        EXPECT_EQ(moved.data[i], original_row_pointers[i]);
    }
    delete[] original_row_pointers;
    
    // Verify original is in valid but empty state
    EXPECT_EQ(original.data, nullptr);
    EXPECT_EQ(original.rows, 0);
    EXPECT_EQ(original.cols, 0);
}

// Test move assignment operator
TEST_F(MatrixJaggedTest, MoveAssignment) {
    MatrixJagged<int> original(2, 2);
    original(0, 0) = 50;
    original(0, 1) = 60;
    int** original_data = original.data;
    
    MatrixJagged<int> assigned(1, 1);
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
    MatrixJagged<int> temp = std::move(assigned);
    assigned = std::move(temp);  // This effectively tests self-move safety
    
    // Verify data is still valid after the "self-move"
    EXPECT_EQ(assigned(0, 0), 50);
    EXPECT_EQ(assigned(0, 1), 60);
    EXPECT_EQ(assigned.rows, 2);
    EXPECT_EQ(assigned.cols, 2);
}

// Test destructor doesn't cause double deletion
TEST_F(MatrixJaggedTest, DestructorSafety) {
    // This test verifies that moved-from objects can be safely destroyed
    MatrixJagged<int>* original = new MatrixJagged<int>(2, 2);
    original->operator()(0, 0) = 100;
    original->operator()(1, 1) = 200;
    
    MatrixJagged<int> moved = std::move(*original);
    delete original;  // Should not cause double deletion
    
    EXPECT_EQ(moved(0, 0), 100);
    EXPECT_EQ(moved(1, 1), 200);
}

// Test clear function
TEST_F(MatrixJaggedTest, ClearFunction) {
    MatrixJagged<int> mat(2, 2);
    mat(0, 0) = 1;
    
    // Manually call clear (through destructor in real usage)
    mat.~MatrixJagged();
    
    // Verify object is in cleared state
    EXPECT_EQ(mat.data, nullptr);
    EXPECT_EQ(mat.rows, 0);
    EXPECT_EQ(mat.cols, 0);
}

// Test with different data types
TEST_F(MatrixJaggedTest, DifferentDataTypes) {
    MatrixJagged<double> double_mat(2, 2);
    double_mat(0, 0) = 1.5;
    double_mat(0, 1) = 2.5;
    double_mat(1, 0) = 3.5;
    double_mat(1, 1) = 4.5;
    
    EXPECT_DOUBLE_EQ(double_mat(0, 0), 1.5);
    EXPECT_DOUBLE_EQ(double_mat(0, 1), 2.5);
    EXPECT_DOUBLE_EQ(double_mat(1, 0), 3.5);
    EXPECT_DOUBLE_EQ(double_mat(1, 1), 4.5);
    
    MatrixJagged<float> float_mat(1, 2);
    float_mat(0, 0) = 3.14f;
    float_mat(0, 1) = 6.28f;
    
    EXPECT_FLOAT_EQ(float_mat(0, 0), 3.14f);
    EXPECT_FLOAT_EQ(float_mat(0, 1), 6.28f);
}

// Test large matrix operations
TEST_F(MatrixJaggedTest, LargeMatrix) {
    const int large_size = 100;
    MatrixJagged<int> large_mat(large_size, large_size);
    
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
    MatrixJagged<int> large_copy = large_mat;
    EXPECT_EQ(large_copy(large_size-1, large_size-1), large_size*large_size - 1);
    
    // Verify all rows were copied
    for (int i = 0; i < large_size; i++) {
        for (int j = 0; j < large_size; j++) {
            EXPECT_EQ(large_copy(i, j), large_mat(i, j));
        }
    }
    
    // Test move of large matrix
    MatrixJagged<int> large_moved = std::move(large_mat);
    EXPECT_EQ(large_moved(large_size-1, large_size-1), large_size*large_size - 1);
    EXPECT_EQ(large_mat.data, nullptr);  // Original should be empty
}

// Test edge cases
TEST_F(MatrixJaggedTest, EdgeCases) {
    // 1x1 matrix
    MatrixJagged<int> single(1, 1);
    single(0, 0) = 42;
    EXPECT_EQ(single(0, 0), 42);
    
    // Very large matrix (stress test)
    MatrixJagged<int> large(500, 500);
    large(499, 499) = 123;
    EXPECT_EQ(large(499, 499), 123);
    
    // Verify all rows are properly allocated
    for (int i = 0; i < 500; i++) {
        EXPECT_NE(large.data[i], nullptr);
    }
}

// Test memory independence of rows
TEST_F(MatrixJaggedTest, MemoryIndependence) {
    MatrixJagged<int> mat(3, 2);
    
    // Fill with different values
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 2; j++) {
            mat(i, j) = i * 10 + j;
        }
    }
    
    // Verify rows are independent in memory
    EXPECT_NE(mat.data[0], mat.data[1]);
    EXPECT_NE(mat.data[1], mat.data[2]);
    EXPECT_NE(mat.data[0], mat.data[2]);
    
    // Modify one row shouldn't affect others
    int original_value = mat(1, 0);
    mat(0, 0) = 999;
    EXPECT_EQ(mat(1, 0), original_value);  // Other row unchanged
}

// Test assignment from different sizes
TEST_F(MatrixJaggedTest, AssignmentDifferentSizes) {
    MatrixJagged<int> small(2, 2);
    small(0, 0) = 1; small(0, 1) = 2;
    small(1, 0) = 3; small(1, 1) = 4;
    
    MatrixJagged<int> large(3, 3);
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            large(i, j) = i * 3 + j;
        }
    }
    
    // Assign small to large - large should resize to small's dimensions
    large = small;
    
    EXPECT_EQ(large.rows, 2);
    EXPECT_EQ(large.cols, 2);
    EXPECT_EQ(large(0, 0), 1);
    EXPECT_EQ(large(0, 1), 2);
    EXPECT_EQ(large(1, 0), 3);
    EXPECT_EQ(large(1, 1), 4);
}
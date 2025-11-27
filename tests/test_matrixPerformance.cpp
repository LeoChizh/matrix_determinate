#include <gtest/gtest.h>
#include "MatrixPerformanceTest.hpp"

TEST(MatrixPerformance, ComprehensiveComparison) {
    // Test with reasonable sizes for meaningful performance data
    const size_t maxSize = 200;
    
    std::cout << "\n\n";
    std::cout << "========================================\n";
    std::cout << "    MATRIX PERFORMANCE COMPARISON\n";
    std::cout << "========================================\n";
    
    MatrixPerformanceTest::runComprehensiveTests<Matrix<int>>(maxSize);
    MatrixPerformanceTest::runComprehensiveTests<Matrix<double>>(maxSize);
    MatrixPerformanceTest::runComprehensiveTests<MatrixJagged<int>>(maxSize);
    MatrixPerformanceTest::runComprehensiveTests<MatrixJagged<double>>(maxSize);
}

TEST(MatrixPerformance, SpecificOperations) {
    // Test specific operations that are most relevant
    std::cout << "\n\n";
    std::cout << "========================================\n";
    std::cout << "    SPECIFIC OPERATIONS COMPARISON\n";
    std::cout << "========================================\n";
    
    // Just test 100x100 for quick comparison
    MatrixPerformanceTest::runComprehensiveTests<Matrix<double>>(100);
    MatrixPerformanceTest::runComprehensiveTests<MatrixJagged<double>>(100);
}
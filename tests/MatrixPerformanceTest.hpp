#pragma once
#include <chrono>
#include <vector>
#include <string>
#include <iostream>
#include <iomanip>
#include <random>
#include <functional>
#include <map>
#include <algorithm>
#include "Matrix.hpp"
#include "MatrixJagged.hpp"
#include "MatrixAlgorithms.hpp"

class MatrixPerformanceTest {
public:
    using Clock = std::chrono::high_resolution_clock;
    using Milliseconds = std::chrono::milliseconds;
    using Microseconds = std::chrono::microseconds;
    
    struct TestResult {
        std::string testName;
        std::string matrixType;
        size_t matrixSize;
        long long durationMicros;
        bool success;
        std::string additionalInfo;
    };

    template<typename MatrixType>
    static void runComprehensiveTests(size_t maxSize = 1000) {
        std::cout << "=== Performance Test for " << getMatrixTypeName<MatrixType>() << " ===\n";
        
        std::vector<TestResult> results;
        
        // Test different sizes - adjusted for reasonable performance
        std::vector<size_t> sizes = {10, 50, 100, 200};
        if (maxSize >= 500) sizes.push_back(500);
        
        for (size_t size : sizes) {
            if (size > maxSize) continue;
            
            std::cout << "\n--- Testing size " << size << "x" << size << " ---\n";
            
            auto result1 = testIdentityCreation<MatrixType>(size);
            results.push_back(result1);
            printSingleResult(result1);
            
            auto result2 = testDeterminantCalculation<MatrixType>(size);
            results.push_back(result2);
            printSingleResult(result2);
            
            auto result3 = testMatrixMultiplication<MatrixType>(size);
            results.push_back(result3);
            printSingleResult(result3);
            
            auto result4 = testMatrixAddition<MatrixType>(size);
            results.push_back(result4);
            printSingleResult(result4);
        }
        
        printSummary(results);
    }

private:
    template<typename MatrixType>
    static std::string getMatrixTypeName() {
        if constexpr (std::is_same_v<MatrixType, Matrix<int>>) return "Matrix<int>";
        if constexpr (std::is_same_v<MatrixType, Matrix<double>>) return "Matrix<double>";
        if constexpr (std::is_same_v<MatrixType, MatrixJagged<int>>) return "MatrixJagged<int>";
        if constexpr (std::is_same_v<MatrixType, MatrixJagged<double>>) return "MatrixJagged<double>";
        return "Unknown";
    }

    // Generic error string function that works with any MatrixType
    template<typename MatrixType>
    static std::string getErrorString(typename MatrixAlgorithms<MatrixType>::DeterminantError error) {
        switch (error) {
            case MatrixAlgorithms<MatrixType>::DeterminantError::SUCCESS: return "SUCCESS";
            case MatrixAlgorithms<MatrixType>::DeterminantError::NOT_SQUARE: return "NOT_SQUARE";
            case MatrixAlgorithms<MatrixType>::DeterminantError::OVERFLOW_NUMBER: return "OVERFLOW";
            default: return "UNKNOWN";
        }
    }

    template<typename MatrixType>
    static TestResult testIdentityCreation(size_t size) {
        TestResult result;
        result.testName = "Identity Creation";
        result.matrixType = getMatrixTypeName<MatrixType>();
        result.matrixSize = size;
        
        try {
            MatrixType reference(size, size);
            
            auto start = Clock::now();
            auto identity = MatrixAlgorithms<MatrixType>::createIdentity(reference);
            auto end = Clock::now();
            
            result.durationMicros = std::chrono::duration_cast<Microseconds>(end - start).count();
            result.success = identity.has_value();
            result.additionalInfo = result.success ? "Success" : "Not square (expected)";
        } catch (const std::exception& e) {
            result.durationMicros = 0;
            result.success = false;
            result.additionalInfo = std::string("Exception: ") + e.what();
        }
        
        return result;
    }

    template<typename MatrixType>
    static TestResult testDeterminantCalculation(size_t size) {
        TestResult result;
        result.testName = "Determinant Calculation";
        result.matrixType = getMatrixTypeName<MatrixType>();
        result.matrixSize = size;
        
        try {
            MatrixType matrix(size, size);
            
            // Use different filling strategies for integer vs floating-point
            if constexpr (std::is_integral_v<typename MatrixType::value_type>) {
                fillVerySmallValues(matrix); // Use very small values for integers
            } else {
                fillWellConditioned(matrix); // Use normal values for floating-point
            }
            
            auto start = Clock::now();
            auto detResult = MatrixAlgorithms<MatrixType>::determinantLU(matrix);
            auto end = Clock::now();
            
            result.durationMicros = std::chrono::duration_cast<Microseconds>(end - start).count();
            result.success = detResult.success();
            
            if (detResult.success()) {
                // Format the determinant value nicely
                if constexpr (std::is_integral_v<typename MatrixType::value_type>) {
                    result.additionalInfo = "Determinant: " + std::to_string(detResult.value);
                } else {
                    char buffer[64];
                    snprintf(buffer, sizeof(buffer), "Determinant: %.2e", detResult.value);
                    result.additionalInfo = buffer;
                }
            } else {
                result.additionalInfo = "Error: " + getErrorString<MatrixType>(detResult.error);
            }
        } catch (const std::exception& e) {
            result.durationMicros = 0;
            result.success = false;
            result.additionalInfo = std::string("Exception: ") + e.what();
        }
        
        return result;
    }

    template<typename MatrixType>
    static TestResult testMatrixMultiplication(size_t size) {
        TestResult result;
        result.testName = "Matrix Multiplication";
        result.matrixType = getMatrixTypeName<MatrixType>();
        result.matrixSize = size;
        
        try {
            MatrixType A(size, size);
            MatrixType B(size, size);
            fillSmallValues(A);
            fillSmallValues(B);
            
            auto start = Clock::now();
            MatrixType C = multiplyMatrices(A, B);
            auto end = Clock::now();
            
            result.durationMicros = std::chrono::duration_cast<Microseconds>(end - start).count();
            result.success = true;
            result.additionalInfo = "Multiplication completed";
        } catch (const std::exception& e) {
            result.durationMicros = 0;
            result.success = false;
            result.additionalInfo = std::string("Exception: ") + e.what();
        }
        
        return result;
    }

    template<typename MatrixType>
    static TestResult testMatrixAddition(size_t size) {
        TestResult result;
        result.testName = "Matrix Addition";
        result.matrixType = getMatrixTypeName<MatrixType>();
        result.matrixSize = size;
        
        try {
            MatrixType A(size, size);
            MatrixType B(size, size);
            fillSmallValues(A);
            fillSmallValues(B);
            
            auto start = Clock::now();
            MatrixType C = addMatrices(A, B);
            auto end = Clock::now();
            
            result.durationMicros = std::chrono::duration_cast<Microseconds>(end - start).count();
            result.success = true;
            result.additionalInfo = "Addition completed";
        } catch (const std::exception& e) {
            result.durationMicros = 0;
            result.success = false;
            result.additionalInfo = std::string("Exception: ") + e.what();
        }
        
        return result;
    }

    template<typename MatrixType>
    static void fillVerySmallValues(MatrixType& matrix) {
        std::random_device rd;
        std::mt19937 gen(rd());
        
        using ValueType = typename MatrixType::value_type;
        
        if constexpr (std::is_integral_v<ValueType>) {
            // Use VERY small values for integers to absolutely avoid overflow
            std::uniform_int_distribution<ValueType> dist(-1, 1);
            for (size_t i = 0; i < matrix.rows; ++i) {
                for (size_t j = 0; j < matrix.cols; ++j) {
                    matrix(i, j) = dist(gen);
                }
            }
            
            // Make diagonally dominant with small values
            for (size_t i = 0; i < matrix.rows; ++i) {
                matrix(i, i) = ValueType{2};
            }
        } else {
            fillWellConditioned(matrix);
        }
    }

    template<typename MatrixType>
    static void fillWellConditioned(MatrixType& matrix) {
        std::random_device rd;
        std::mt19937 gen(rd());
        
        using ValueType = typename MatrixType::value_type;
        
        if constexpr (std::is_integral_v<ValueType>) {
            std::uniform_int_distribution<ValueType> dist(-3, 3);
            for (size_t i = 0; i < matrix.rows; ++i) {
                for (size_t j = 0; j < matrix.cols; ++j) {
                    matrix(i, j) = dist(gen);
                }
            }
            
            // Make diagonally dominant to ensure non-singular
            for (size_t i = 0; i < matrix.rows; ++i) {
                matrix(i, i) = ValueType{10};
            }
        } else {
            std::uniform_real_distribution<ValueType> dist(-5.0, 5.0);
            for (size_t i = 0; i < matrix.rows; ++i) {
                for (size_t j = 0; j < matrix.cols; ++j) {
                    matrix(i, j) = dist(gen);
                }
            }
            
            // Make diagonally dominant
            for (size_t i = 0; i < matrix.rows; ++i) {
                matrix(i, i) = ValueType{20};
            }
        }
    }

    template<typename MatrixType>
    static void fillSmallValues(MatrixType& matrix) {
        std::random_device rd;
        std::mt19937 gen(rd());
        
        using ValueType = typename MatrixType::value_type;
        
        if constexpr (std::is_integral_v<ValueType>) {
            std::uniform_int_distribution<ValueType> dist(-2, 2);
            for (size_t i = 0; i < matrix.rows; ++i) {
                for (size_t j = 0; j < matrix.cols; ++j) {
                    matrix(i, j) = dist(gen);
                }
            }
        } else {
            std::uniform_real_distribution<ValueType> dist(-5.0, 5.0);
            for (size_t i = 0; i < matrix.rows; ++i) {
                for (size_t j = 0; j < matrix.cols; ++j) {
                    matrix(i, j) = dist(gen);
                }
            }
        }
    }

    template<typename MatrixType>
    static MatrixType multiplyMatrices(const MatrixType& A, const MatrixType& B) {
        MatrixType result(A.rows, B.cols);
        
        for (size_t i = 0; i < A.rows; ++i) {
            for (size_t j = 0; j < B.cols; ++j) {
                typename MatrixType::value_type sum{};
                for (size_t k = 0; k < A.cols; ++k) {
                    sum += A(i, k) * B(k, j);
                }
                result(i, j) = sum;
            }
        }
        
        return result;
    }

    template<typename MatrixType>
    static MatrixType addMatrices(const MatrixType& A, const MatrixType& B) {
        MatrixType result(A.rows, A.cols);
        
        for (size_t i = 0; i < A.rows; ++i) {
            for (size_t j = 0; j < A.cols; ++j) {
                result(i, j) = A(i, j) + B(i, j);
            }
        }
        
        return result;
    }

    static void printSingleResult(const TestResult& result) {
        // Use "us" instead of μs to avoid encoding issues
        std::cout << std::left << std::setw(25) << result.testName
                  << std::setw(8) << result.matrixSize << "x" << result.matrixSize
                  << std::setw(10) << result.durationMicros << "us"
                  << std::setw(12) << (result.success ? "PASS" : "FAIL")
                  << result.additionalInfo << "\n";
    }

    static void printSummary(const std::vector<TestResult>& results) {
        std::cout << "\n" << std::string(70, '=') << "\n";
        std::cout << "PERFORMANCE SUMMARY\n";
        std::cout << std::string(70, '=') << "\n";
        
        // Group by matrix type
        std::map<std::string, std::vector<TestResult>> groupedResults;
        for (const auto& result : results) {
            groupedResults[result.matrixType].push_back(result);
        }
        
        for (const auto& [matrixType, typeResults] : groupedResults) {
            std::cout << "\n" << matrixType << ":\n";
            std::cout << std::string(50, '-') << "\n";
            
            for (const auto& result : typeResults) {
                if (result.success) {
                    std::cout << std::left << std::setw(25) << result.testName
                              << std::setw(15) << (std::to_string(result.matrixSize) + "x" + std::to_string(result.matrixSize))
                              << std::setw(12) << (std::to_string(result.durationMicros) + " us")
                              << "\n";
                }
            }
        }
        
        // Calculate and show speed comparisons
        std::cout << "\n" << std::string(70, '=') << "\n";
        std::cout << "SPEED COMPARISON (Lower is better)\n";
        std::cout << std::string(70, '=') << "\n";
        
        // Compare all matrix types for common operations
        compareAllMatrixTypes(results, "Matrix Multiplication", 100);
        compareAllMatrixTypes(results, "Matrix Addition", 100);
        compareAllMatrixTypes(results, "Determinant Calculation", 100);
    }

    static void compareAllMatrixTypes(const std::vector<TestResult>& results, 
                                    const std::string& testName, size_t size) {
        std::cout << "\n" << testName << " (" << size << "x" << size << "):\n";
        
        std::map<std::string, long long> times;
        for (const auto& result : results) {
            if (result.testName == testName && result.matrixSize == size && result.success) {
                times[result.matrixType] = result.durationMicros;
            }
        }
        
        if (times.size() > 1) {
            // Find the fastest
            auto fastest = std::min_element(times.begin(), times.end(),
                [](const auto& a, const auto& b) { return a.second < b.second; });
            
            for (const auto& [matrixType, time] : times) {
                double ratio = static_cast<double>(time) / fastest->second;
                std::cout << "  " << std::setw(20) << matrixType 
                          << std::setw(8) << time << " us"
                          << " (" << std::fixed << std::setprecision(2) << ratio << "x)\n";
            }
        } else if (times.size() == 1) {
            std::cout << "  " << std::setw(20) << times.begin()->first 
                      << std::setw(8) << times.begin()->second << " us (only one type succeeded)\n";
        } else {
            std::cout << "  No successful results for this test\n";
        }
    }
};
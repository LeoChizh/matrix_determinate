cmake_minimum_required(VERSION 3.15)
project(SearchingTree VERSION 1.0.0 LANGUAGES CXX)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Set compiler warnings
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wpedantic -Werror")
endif()

# Create a library for your core functionality

#add_library(ExampleLib)

# Target include directories

#target_include_directories(ExampleLib PUBLIC include)

# Main executable
add_executable(${PROJECT_NAME}
    src/main.cpp
)

# Link main executable with the library

#target_link_libraries(${PROJECT_NAME} ExampleLib)

# Add tests subdirectory if it exists
if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/tests)
    enable_testing()
    add_subdirectory(tests)
endif()
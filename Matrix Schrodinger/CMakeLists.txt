cmake_minimum_required(VERSION 3.23)
project(Matrix_Schrodinger)

set(CMAKE_CXX_STANDARD 17)

add_subdirectory(eigen)
set(EIGEN3_INCLUDE_DIR "${CMAKE_CURRENT_LIST_DIR}/eigen/Eigen")
include_directories(${EIGEN3_INCLUDE_DIR})

add_executable(Matrix_Schrodinger Source.cpp)
target_link_libraries (Matrix_Schrodinger Eigen3::Eigen)

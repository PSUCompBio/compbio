cmake_minimum_required(VERSION 2.6)
project(compbio)
enable_language(CXX)

# Eigen support
option (USE_EIGEN  "Build library with Eigen support" ON)
set (EIGEN_INCLUDE_DIR "INCLUDE directory for Eigen library" third-party-libs/eigen3)
include_directories(${EIGEN_INCLUDE_DIR})

#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fpic")

# All the source code we care about
add_subdirectory(src)


cmake_minimum_required(VERSION 3.10)
project(BVPSolver)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

find_program(GNUPLOT_EXECUTABLE gnuplot REQUIRED)

add_executable(bvp_demo main.cpp)

# Add gnuplot path to compile definitions
add_compile_definitions(GNUPLOT_PATH="${GNUPLOT_EXECUTABLE}")

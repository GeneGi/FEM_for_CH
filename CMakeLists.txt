cmake_minimum_required(VERSION 3.7)
project(code)

set(CMAKE_CXX_STANDARD 11)

set(SOURCE_FILES src/generate_info_matrix.cpp
        src/gauss_quadrature.cpp
        src/basis_function.cpp
        src/assemble_matrix.cpp
        src/simple_function.cpp
        src/load_vector.cpp
        src/fe_solution_triangle.cpp
        src/initial_function.cpp
        src/treat_boundary.cpp
        src/boundary_function.cpp
        src/write_file.cpp)

add_library(fem ${SOURCE_FILES})

add_executable(code src/main.cpp)
target_link_libraries(code fem)

add_executable(test_solve test/test_solve.cpp)
target_link_libraries(test_solve fem)

add_executable(test_basis_function test/test_basis_function.cpp)
target_link_libraries(test_basis_function fem)

add_executable(test_gauss test/test_gauss_quadrature.cpp)
target_link_libraries(test_gauss fem)

add_executable(read_from_file test/read_from_file.cpp)
target_link_libraries(read_from_file fem)

add_executable(test_openmp test/test_openmp.cpp)
cmake_minimum_required(VERSION 3.7)

project(code)

set(CMAKE_CXX_COMPILER "/usr/local/bin/g++-6")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O0 -fopenmp")

LIST(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake") # Add the directory where FindSuiteSparse.cmake module can be found.

set(SuiteSparse_USE_LAPACK_BLAS ON)
find_package(SuiteSparse QUIET NO_MODULE)  # 1st: Try to locate the *config.cmake file.
if(NOT SuiteSparse_FOUND)
    set(SuiteSparse_VERBOSE ON)
    find_package(SuiteSparse REQUIRED) # 2nd: Use FindSuiteSparse.cmake module
    include_directories(${SuiteSparse_INCLUDE_DIRS})
else()
    message(STATUS "Find SuiteSparse : include(${USE_SuiteSparse})")
    include(${USE_SuiteSparse})
endif()
MESSAGE(STATUS "SuiteSparse_LIBS: ${SuiteSparse_LIBRARIES}")

#find_package(OpenMP)
#if(OPENMP_FOUND)
#    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
#    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
#    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
#endif()


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
        src/write_file.cpp
        src/read_mesh_from_file.cpp)

add_library(fem ${SOURCE_FILES})

add_executable(code src/main.cpp)
target_link_libraries(code fem)

add_executable(test_solve test/test_solve.cpp)
target_link_libraries(test_solve fem ${SuiteSparse_LIBRARIES})

add_executable(test_basis_function test/test_basis_function.cpp)
target_link_libraries(test_basis_function fem)

add_executable(test_gauss test/test_gauss_quadrature.cpp)
target_link_libraries(test_gauss fem)

add_executable(read_from_file test/read_from_file.cpp)
target_link_libraries(read_from_file fem)

add_executable(test_openmp test/test_openmp.cpp)

add_executable(test_speed test/test_speed.cpp)
target_link_libraries(test_speed fem)

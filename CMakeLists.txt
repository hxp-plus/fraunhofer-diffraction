cmake_minimum_required(VERSION 3.17)
project(fraunhofer-diffraction)

set(CMAKE_C_STANDARD 99)
add_executable(fraunhofer-diffraction main.c)

target_link_libraries(fraunhofer-diffraction gsl gslcblas m)
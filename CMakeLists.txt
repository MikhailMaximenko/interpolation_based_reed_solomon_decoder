# CMakeList.txt : Top-level CMake project file, do global configuration
# and include sub-projects here.
#
cmake_minimum_required (VERSION 3.8)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

project ("fast_interpolation_reed_solomon_decoder")

# Include sub-projects.
add_subdirectory("fast_interpolation_reed_solomon_decoder")

add_subdirectory("tests")
add_executable(decoder main.cpp)
target_link_libraries(decoder PUBLIC fast_interpolation_reed_solomon_decoder)

target_include_directories(decoder PUBLIC "${PROJECT_BINARY_DIR}" "${PROJECT_SOURCE_DIR}/fast_interpolation_reed_solomon_decoder")
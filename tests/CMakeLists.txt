cmake_minimum_required(VERSION 3.21)
project(decoder-test)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(TSOURCES galois-field-test.cpp)
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
include(FetchContent)
FetchContent_Declare(
    googletest
    GIT_REPOSITORY https://github.com/google/googletest.git
    GIT_TAG        v1.15.2
)
FetchContent_MakeAvailable(googletest)

add_executable(tests ${TSOURCES})

target_link_libraries(tests
    PUBLIC
        fast_interpolation_reed_solomon_decoder
        GTest::gtest_main
)

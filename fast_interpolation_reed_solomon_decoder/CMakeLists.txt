file(GLOB HEADERS "*.h")
file(GLOB SOURCES "*.cpp")

add_library(fast_interpolation_reed_solomon_decoder ${HEADERS} ${SOURCES})
target_include_directories(fast_interpolation_reed_solomon_decoder PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})


# CMakeList.txt : CMake project for fast_interpolation_reed_solomon_decoder, include source and define
# project specific logic here.
#

# Add source to this project's executable.
add_executable (fast_interpolation_reed_solomon_decoder 
	"fast_interpolation_reed_solomon_decoder.cpp" "fast_interpolation_reed_solomon_decoder.h" 
	"galois_field.h" "galois_field.cpp" 
	"fft.h" "FFT7.CPP" "FFT63.CPP" "FFT127.CPP" "FFT255.CPP" "FFT511.CPP")

if (CMAKE_VERSION VERSION_GREATER 3.12)
  set_property(TARGET fast_interpolation_reed_solomon_decoder PROPERTY CXX_STANDARD 20)
endif()

# TODO: Add tests and install targets if needed.

#pragma once

#include <vector>

void FFT7(unsigned*, unsigned*);
void FFT63(unsigned*, unsigned*);
void FFT127(unsigned*, unsigned*);
void FFT255(unsigned*, unsigned*);
void FFT511(unsigned*, unsigned*);


void inline call_fft(std::vector<unsigned>& src, std::vector<unsigned>& dst, unsigned length, size_t& additions, size_t& multiplications) {
	if (length == 7) {
		FFT7(src.data(), dst.data());
		multiplications += 6;
		additions += 25;
		return;
	}
	if (length == 63) {
		FFT63(src.data(), dst.data());
		multiplications += 97;
		additions += 872;
		return;
	}
	if (length == 127) {
		FFT127(src.data(), dst.data());
		multiplications += 216;
		additions += 2839;
		return;
	}
	if (length == 255) {
		FFT255(src.data(), dst.data());

		return;
	}
	if (length == 511) {
		FFT511(src.data(), dst.data());
		return;
	}
}
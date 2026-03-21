#pragma once

#include <vector>

void FFT7(unsigned*, unsigned*);
void FFT63(unsigned*, unsigned*);
void FFT127(unsigned*, unsigned*);
void FFT255(unsigned*, unsigned*);
void FFT511(unsigned*, unsigned*);


void inline call_fft(std::vector<unsigned>& src, std::vector<unsigned>& dst) {
	if (src.size() <=64) {
		FFT63(src.data(), dst.data());
		return;
	}
	//if (src.size() <= 63) {
	//	FFT63(src.data(), dst.data());
	//	return;
	//}
	//if (src.size() <= 127) {
	//	FFT127(src.data(), dst.data());
	//	return;
	//}
	//if (src.size() <= 255) {
	//	FFT255(src.data(), dst.data());
	//	return;
	//}
	//if (src.size() <= 511) {
	//	FFT511(src.data(), dst.data());
	//	return;
	//}
}
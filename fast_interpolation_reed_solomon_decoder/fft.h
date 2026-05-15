#pragma once

#include <vector>
#include<iostream>
#include "galois_field.h"
void FFT7(unsigned*, unsigned*);
void FFT63(unsigned*, unsigned*);
void FFT127(unsigned*, unsigned*);
void FFT255(unsigned*, unsigned*);
void FFT511(unsigned*, unsigned*);
void FFT_1023(const unsigned*, unsigned*, const galois_field&);
void FFT_2047(const unsigned*, unsigned*, const galois_field&);
void FFT_4095(const unsigned*, unsigned*, const galois_field&);


void inline call_fft(galois_field &gf, std::vector<unsigned>& src, std::vector<unsigned>& dst, unsigned length) {
		std::cout << "here\n";
	if (length == 7) {
		FFT7(src.data(), dst.data());
		gf._multiplications += 6;
		gf._additions += 25;
		return;
	}
	if (length == 63) {
		FFT63(src.data(), dst.data());
		gf._multiplications += 97;
		gf._additions += 872;
		return;
	}
	if (length == 127) {
		FFT127(src.data(), dst.data());
		gf._multiplications += 216;
		gf._additions += 2839;
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
	if (length == 1023) {
		FFT_1023(src.data(), dst.data(), gf);
	}
	if (length == 2047) {
		FFT_2047(src.data(), dst.data(), gf);
	}
	if (length == 4095) {
		FFT_4095(src.data(), dst.data(), gf);
	}
}
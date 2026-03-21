// fast_interpolation_reed_solomon_decoder.h : Include file for standard system include files,
// or project specific include files.

#pragma once
#include "galois_field.h"
#include <iostream>
#include <array>
#include<random>
struct InterpolationBasedFastRSDecoder {
	galois_field _gf;

	size_t _n, _k, _t;
	std::array<std::vector<unsigned>, 12> _tmp;

	InterpolationBasedFastRSDecoder(galois_field const&, unsigned, unsigned);

	std::vector<unsigned> encode(std::vector<unsigned> &);

	void decode(std::vector<unsigned>&);

};

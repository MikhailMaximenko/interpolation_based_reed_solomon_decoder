// fast_interpolation_reed_solomon_decoder.h : Include file for standard system include files,
// or project specific include files.

#pragma once
#include "galois_field.h"
#include <iostream>
#include <array>

struct InterpolationBasedFastRSDecoder {
	galois_field _gf;

	size_t _n, _k, _t;
	std::array<std::vector<size_t>, 12> _tmp;

	std::vector<size_t> encode(std::vector<size_t> const&);

	void decode(std::vector<size_t>&);

};

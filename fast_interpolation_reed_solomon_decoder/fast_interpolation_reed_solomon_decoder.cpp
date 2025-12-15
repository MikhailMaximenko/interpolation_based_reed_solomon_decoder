// fast_interpolation_reed_solomon_decoder.cpp : Defines the entry point for the application.
//

#include "fast_interpolation_reed_solomon_decoder.h"
#include "galois_field.h"
#include <bit>



void InterpolationBasedFastRSDecoder::decode(std::vector<size_t>& cw) {
	_gf.IDFT(cw, _tmp[0], _gf._exp_table[1]);
	_tmp[1].insert(_tmp[1].begin(), _tmp[0].rbegin(), _tmp[0].rend() + (_n - 2 * _t));
	_tmp[2].insert(_tmp[2].begin(), _tmp[0].rbegin() - _t, _tmp[0].rend() + (_n - 2 * _t - 1)); // n - 2 * t + 1 ..n - t
	_gf.SOLVE_TOEPITZ(_tmp[1], _tmp[2], _t, _tmp[3]);



	using std::swap;
	swap(cw, _tmp[1]);
	_gf.sub_poly(_tmp[1], _gf.DFT(), cw);


}

int main()
{
	galois_field gf(11, 1, { 9, 1 });
	for (auto i : gf._log_table) {
		std::cout << i << " ";
	}
	std::cout << "\n";
	for (auto i : gf._exp_table) {
		std::cout << i << " ";
	}
	std::cout << "\n";

	std::cout << gf.divide(8, 10) << "\n";


	return 0;
}

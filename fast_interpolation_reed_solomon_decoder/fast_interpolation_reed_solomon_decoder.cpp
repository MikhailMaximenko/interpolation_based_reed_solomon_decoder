// fast_interpolation_reed_solomon_decoder.cpp : Defines the entry point for the application.
//

#include "fast_interpolation_reed_solomon_decoder.h"
#include "galois_field.h"
#include "fft.h"
#include <bit>



void InterpolationBasedFastRSDecoder::decode(std::vector<unsigned>& cw) {
	_gf.IDFT(cw, _tmp[0]);
	_tmp[1].insert(_tmp[1].begin(), _tmp[0].rbegin(), _tmp[0].rend() + (_n - 2 * _t));
	_tmp[2].insert(_tmp[2].begin(), _tmp[0].rbegin() - _t, _tmp[0].rend() + (_n - 2 * _t - 1)); // n - 2 * t + 1 ..n - t
	_gf.SOLVE_TOEPITZ(_tmp[1], _tmp[2], _t, _tmp[3]);
	for (ptrdiff_t i = _k - 1; i >= 0; --i) {
		_tmp[0][i] = 0;
		for (size_t j = 0; j < _t; ++j) {
			_tmp[0][i] = _gf.add(_tmp[0][i], _gf.multiply(_tmp[0][i + j], _tmp[3][j]));
		}
	}

	using std::swap;
	swap(cw, _tmp[1]);
	_gf.sub_poly(_tmp[1], _gf.DFT(_tmp[0], _tmp[2]), cw);


}

int main()
{
	galois_field gf2(3, 0xb, 3);
	std::cout << "here\n";
	for (size_t i = 0; i < (1 << 3); ++i) {
		std::cout << gf2._exp_table[i] << " ";
	}
	std::cout << "\n";
	for (size_t i = 0; i < (1 << 3); ++i) {
		std::cout << gf2._log_table[i] << " ";
	}

	std::cout << "\n";

	std::vector<unsigned> a{ 2, 3, 0, 0, 0, 0, 0 };
	std::vector<unsigned> b{ 1, 1, 1, 1, 0, 0, 0 };
	std::vector<unsigned> c(7);
	gf2.fast_poly_multiplication(a, b, c);
	//for (auto i : gf2._a_tmp) {
	//	std::cout << i << " ";
	//}
	//std::cout << "\n";
	//for (auto i : gf2._b_tmp) {
	//	std::cout << i << " ";
	//}
	//std::cout << "\n";
	//call_fft(gf2._b_tmp, gf2._a_tmp);
	//for (auto i : gf2._a_tmp) {
		//std::cout << i << " ";
	//}
	//std::cout << "\n";
	for (auto i : c) {
		std::cout << i << " ";
	}
	std::cout << "\n";

	std::vector<unsigned> aa = { 1, 2, 3, 4, 5, 6, 7 };
	std::vector<unsigned> ab(7);
	call_fft(aa, ab);
	aa = ab;
	call_fft(aa, ab);
	for (auto i : ab) {
		std::cout << i << " ";
	}
	std::cout << "\n";

	std::cout << "snd\n";
	std::vector<unsigned> a1 = { 2, 3, 0, 0, 0, 0, 0, 0 }, b1 = {2, 3, 0, 0, 0, 0, 0, 0}, f, g;
	f.resize(8); g.resize(8);
	gf2.DFT(a1, f, 4, 0);
	std::cout << "fst\n";
	gf2.DFT(b1, g, 4, 0);

	std::cout << "completed!\n";

	for (size_t i = 0; i < a1.size(); ++i) {
		a1[i] = gf2.multiply(f[i], g[i]);
		std::cout << a1[i] << " ";
	}
	std::cout << "\n";

	gf2.IDFT(a1, b1, 4, 0);
	for (auto i : b1) {
		std::cout << i << " ";
	}
	std::cout << "\n";
	a1 = { 2, 3, 0, 0, 0, 0, 0, 0 };
	gf2.DFT(a1, f, 4, 0);
	std::cout << "straight\n";
	for (auto i : f) {
		std::cout << i << " ";
	}
	std::cout << "\n";
	std::cout << "inverse\n";
	for (auto i : gf2.IDFT(f, g, 4, 0)) {
		std::cout << i << " ";
	}
	std::cout << "\n";

	std::cout << "id expected for:\n";
	a1 = { 2, 3, 0, 0, 0, 0, 0, 0 };
	for (auto i : a1) {
		std::cout << i << " ";
	}
	std::cout << "\n";
	gf2.DFT(a1, b1, 4, 0);
	gf2.IDFT(b1, a1, 4, 0);
	std::cout << "got:\n";
	for (auto i : a1) {
		std::cout << i << " ";
	}
	std::cout << "\n";

	return 0;
}

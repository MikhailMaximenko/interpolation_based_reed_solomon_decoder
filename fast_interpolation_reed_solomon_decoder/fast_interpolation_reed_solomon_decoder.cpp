// fast_interpolation_reed_solomon_decoder.cpp : Defines the entry point for the application.
//

#include "fast_interpolation_reed_solomon_decoder.h"
#include "galois_field.h"
#include "fft.h"
#include <bit>

InterpolationBasedFastRSDecoder::InterpolationBasedFastRSDecoder(galois_field const& gf, unsigned n, unsigned k)
	: _gf(gf)
	, _n(n)
	, _k(k)
	, _t((n - k) / 2)
{
	for (auto& v : _tmp) {
		v.resize(_n + 1);
	}
}


std::vector<unsigned> InterpolationBasedFastRSDecoder::encode(std::vector<unsigned> & m) {
	std::vector<unsigned> res(_n);
	_gf.DFT(m, res);
	return res;
}


void InterpolationBasedFastRSDecoder::decode(std::vector<unsigned>& cw) {
	if (_t != 0) {
		for (auto& tmp : _tmp) {
			std::fill(tmp.begin(), tmp.end(), 0);
		}
		_gf.IDFT(cw, _tmp[0]);
		std::copy(_tmp[0].begin() + _n - 2 * _t + 1, _tmp[0].begin() + _n, _tmp[1].begin());
		std::copy(_tmp[0].begin() + _n - 2 * _t, _tmp[0].begin() + _n - _t, _tmp[2].begin());
		_gf.SOLVE_TOEPITZ(_tmp[1], _tmp[2], _t - 1, _tmp[3]);

		for (ptrdiff_t i = _k - 1; i >= 0; --i) {
			_tmp[0][i] = 0;
			for (size_t j = 1; j <= _t; ++j) {
				_tmp[0][i] = _gf.add(_tmp[0][i], _gf.multiply(_tmp[0][i + j], _tmp[3][_t - j]));
			}
		}

		using std::swap;
		swap(cw, _tmp[1]);
		_gf.DFT(_tmp[0], _tmp[2]);
		_gf.sub_poly(_tmp[1], _tmp[2], cw);
	}



}

std::random_device rd;
std::mt19937 gen(rd());


std::vector<unsigned> generate_message(unsigned n, unsigned k) {
	std::vector<unsigned> msg(n);
	std::uniform_int_distribution<unsigned> distr(0, n);
	std::generate_n(msg.begin(), k, [&]() { return distr(gen); });
	return msg;
}

std::vector<unsigned> generate_errors(unsigned n, unsigned t) {
	std::vector<unsigned> errors(n);
	std::uniform_int_distribution<unsigned> distr(1, n);
	std::generate_n(errors.begin(), t, [&]() { return distr(gen); });
	std::shuffle(errors.begin(), errors.end(), gen);
	return errors;
}



void test_decoder(galois_field & gf, unsigned n, unsigned k, unsigned iters) {
	InterpolationBasedFastRSDecoder decoder(gf, n, k);
	for (size_t t = 0; t <= (n - k) / 2; ++t) {
		std::cout << "testing t: " << t << "\n";
		decoder._t = t;
		for (size_t _ = 0; _ < iters; ++_) {
			auto msg = generate_message(n, k);
			auto encoded = decoder.encode(msg);
			auto errors = generate_errors(n, t);
			std::vector<unsigned> msg_with_errors(n + 1);
			gf.add_poly(encoded, errors, msg_with_errors, 0);
			decoder.decode(msg_with_errors);
			for (size_t i = 0; i < n; ++i) {
				if (encoded[i] != msg_with_errors[i]) {
					std::cout << "decoding error occured with t:"<< t << "\n";
					gf.print_poly(msg);
					gf.print_poly(encoded);
					gf.print_poly(errors);
					gf.print_poly(msg_with_errors);
					return;
				}
			}
		}
	}
	std::cout << "tests passed\n";


}

int main()
{
	//galois_field gf2(7, 0x83, 7);

	//test_decoder(gf2, 127, 11, 15);

	galois_field gf(3, 0xb, 3);

	std::vector<unsigned> a(20), b(20), c(80), d, e(80);
	a[0] = 1;
	a[8] = 1;
	b[0] = 1;
	b[8] = 1;
	
	d = { 1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4 };

	//gf.add_subpoly_with_modular_shift(d, d, 0, 4, 8, 2, 5);


	//gf.SCHONHAGE_DFT(a, c, 3, 3, 3, 6, 0);
	//gf.SCHONHAGE_DFT(c, e, 3, 6, 3, 6, 0);

	//gf.add_subpoly_with_modular_shift(d, d, 0, 3, 6, 3, );

	gf.SCHONHAGE_STRASSEN_FFT(a, b, c, 9, 0);
	gf.print_poly(c);
	//gf.print_poly(e);
	//gf.print_poly(a);

	return 0;
}

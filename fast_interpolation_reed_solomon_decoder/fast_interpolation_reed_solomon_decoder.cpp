// fast_interpolation_reed_solomon_decoder.cpp : Defines the entry point for the application.
//

#include "fast_interpolation_reed_solomon_decoder.h"
#include "galois_field.h"
#include "fft.h"
#include "berlecamp_massey_decoder.h"
#include <cassert>
#include <bit>

InterpolationBasedFastRSDecoder::InterpolationBasedFastRSDecoder(galois_field const& gf, unsigned n, unsigned k)
	: _gf(gf)
	, _n(n)
	, _k(k)
	, _t((n - k) / 2)
{
	size_t sz = (_n + 1);
	for (auto& v : _tmp) {
		v.resize(sz);
	}
	for (auto& v : _emgcd_tmp) {
		v.first.resize(sz);
		v.second.resize(sz);
	}
}


std::vector<unsigned> InterpolationBasedFastRSDecoder::encode(std::vector<unsigned> & m) {
	std::vector<unsigned> res(_n);
	_gf.DFT(m, res);
	return res;
}


void InterpolationBasedFastRSDecoder::decode(std::vector<unsigned>& cw) {
	// g <=> tmp[0]
	_gf.print_poly(cw);
	_gf.IDFT(cw, _tmp[0]);
	if (_gf.degree(_tmp[0]) <= _k) {
		return ;
	}
	for (auto& tmp : _tmp) {
		std::fill(tmp.begin(), tmp.end(), 0);
	}
	_gf.IDFT(cw, _tmp[0]);

	//_tmp[2][0] = 1;
	//_tmp[2][_n] = 1;
	size_t max_t = (_n - _k) / 2;
	std::copy(_tmp[0].begin() + _n - 2 * max_t + 1, _tmp[0].begin() + _n, _tmp[1].begin());

	//std::copy(_tmp[0].begin() + _n - max_t, _tmp[0].begin() + _n, _tmp[2].begin());
	//std::reverse(_tmp[2].begin(), _tmp[2].begin() + max_t);
	std::reverse(_tmp[1].begin(), _tmp[1].begin() + 2 * max_t);
	//_tmp[1][max_t - 1] = 1;	
	_tmp[2][2 * max_t] = 1;
	std::cout << "here1\n";
	_gf.EMGCD(_tmp[2], _tmp[1], _emgcd_tmp, 0);
	std::cout << "here2\n";

	//_gf.print_poly(_tmp[0]);
	//_gf.print_poly(_tmp[1]);
	//_gf.print_poly(_tmp[2]);
	//_gf.print_poly(_tmp[3]);
	size_t fst_deg = _gf.degree(_emgcd_tmp[0].first);
	size_t t = _gf.degree(_emgcd_tmp[2].second);
	/*if (fst_deg < max_t) {
		t = fst_deg;
	}*/
	std::cout << t << " " << _t << "\n";
	assert(t == _t);
	if (t < _t) {
		_gf.print_poly(_tmp[0]);
		_gf.print_poly(_tmp[1]);
		_gf.print_poly(_emgcd_tmp[0].first);
		_gf.print_poly(_emgcd_tmp[0].second);
		_gf.print_poly(_emgcd_tmp[1].second);
		_gf.print_poly(_emgcd_tmp[2].second);
	}
	std::fill(_tmp[3].begin(), _tmp[3].end(), 0);
	std::fill(_tmp[2].begin(), _tmp[2].end(), 0);
	//std::fill(_tmp[2].begin(), _tmp[2].end(), 0);
	std::fill(_tmp[1].begin(), _tmp[1].end(), 0);

	if (_t != 0) {
		std::copy(_tmp[0].begin() + _n - 2 * _t + 1, _tmp[0].begin() + _n, _tmp[1].begin());
		std::copy(_tmp[0].begin() + _n - 2 * _t, _tmp[0].begin() + _n - _t, _tmp[2].begin());
		_gf.SOLVE_TOEPITZ(_tmp[1], _tmp[2], _t - 1, _tmp[3]);

		// tmp[3] <=> eta
		std::copy(_tmp[3].begin(), _tmp[3].begin() + _t, _tmp[4].begin());
		std::copy(_tmp[0].begin() + _k, _tmp[0].begin() + _k + _t, _tmp[5].begin());
		std::reverse(_tmp[4].begin(), _tmp[4].begin() + _t + 1); // D(x)
		//std::cout << "before reverse:\n";
		//std::cout << "k: " << _k << " t: " << _t << "\n";
		//_gf.print_poly(_tmp[0]);
		//_gf.print_poly(_tmp[5]);
		std::reverse(_tmp[5].begin(), _tmp[5].begin() + _t); // h0 .. h_(t-1)
		_tmp[4][0] = 1;
		_gf.inv_poly(_tmp[4], _tmp[6], _k + _t);
		_gf.remainder_of_power(_gf.fast_poly_multiplication(_tmp[5], _tmp[4], _tmp[8]), _t);
		_gf.remainder_of_power(_gf.fast_poly_multiplication(_tmp[6], _tmp[8], _tmp[7]), _k + _t);
		std::reverse(_tmp[7].begin(), _tmp[7].begin() + _t + _k);
		std::copy(_tmp[7].begin(), _tmp[7].begin() + _k, _tmp[0].begin());
		//std::cout << "got:\n";
			//_gf.print_poly(_tmp[3]);
			//_gf.print_poly(_tmp[4]);
			//_gf.print_poly(_tmp[5]);
			//_gf.print_poly(_tmp[8]);
			//_gf.print_poly(_tmp[6]);
			//_gf.print_poly(_tmp[7]);
		//std::cout << "expected:\n";


		 

		//for (ptrdiff_t i = _k - 1; i >= 0; --i) {
		//	_tmp[0][i] = 0; 
		//	for (size_t j = 1; j <= _t; ++j) {
		//		_tmp[0][i] = _gf.add(_tmp[0][i], _gf.multiply(_tmp[0][i + j], _tmp[3][_t - j]));
		//	}
		//}
		//_gf.print_poly(_tmp[0]);

		using std::swap;
		swap(cw, _tmp[1]);
		_gf.DFT(_tmp[0], _tmp[2]);
		_gf.add_poly(_tmp[1], _tmp[2], cw, 0);
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
	encoding::bch_decoder decoder(gf, n, k, 2, n - k + 1, 1);
	for (size_t t = 0; t <= (n - k)/2; ++t) {
		std::cout << "testing t: " << t << "\n";
		//decoder._t = t;
		for (size_t _ = 0; _ < iters; ++_) {
			auto msg = generate_message(n, k);
			std::vector<unsigned> encoded(n);
			gf.DFT(msg, encoded);
			//auto encoded = decoder.encode(msg);
			auto errors = generate_errors(n, t);
			std::vector<unsigned> msg_with_errors(n);
			gf.add_poly(encoded, errors, msg_with_errors, 0);
			size_t sz = n * gf._m;
			linalg::bit_vector vt(sz, false);
			decoder._gf.reset_counters();
			gf.translate_to_bit_vector(msg_with_errors, vt);
			decoder._decoding.cp(vt);
			auto res = decoder.decode();
			gf.translate_bit_vector(res, msg_with_errors);

			std::cout << "additions: " << decoder._gf._additions << " multiplications: " << decoder._gf._multiplications << "\n";
			decoder._gf.reset_counters();
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
		//break;
	}
	std::cout << "tests passed\n";


}

int main()
{
	//galois_field gf2(7, 0x83, 7);

	//test_decoder(gf2, 127, 11, 15);

	//galois_field gf(3, 0xb, 3);

	// 2 2
	galois_field gf(3, 0xb, 3);

	test_decoder(gf, 7, 1, 10);
	//test_decoder(gf, 511, 256, 10);

	//std::vector<unsigned> a(7), b(7), c(80), d(80), e(80), f(80);
	//a[1] = 1;
	//gf.DFT(a, b);
	//gf.print_poly(b);
	//a = { 38, 45, 31, 40, 3, 55, 51, 6, 50, 30, 31, 62, 5, 44, 38, 42, 18, 33, 47, 30, 39, 11, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
	//a[0] = 2;
	//a[1] = 2;
	//a[2] = 1;
	//a[3] = 7;
	//a[4] = 2;
	//b[0] = 3;
	//b[1] = 1;
	//c[0] = 1;
	//c[1] = 3;
	//c[2] = 1;
	// expected result 1 7
	//gf.inv_poly(c, d, 4);
	//gf.remainder_of_power(gf.fast_poly_multiplication(a, c, f), 2);
	//gf.print_poly(f);
	//gf.fast_poly_multiplication(f, d, e);
	//gf.print_poly(d);
	//gf.print_poly(e);
	////a[0] = 3;
	////a[1] = 63;
	//b[0] = 56;
	//b[2] = 60;
	//gf.remainder_of_power(a, 8);
	//gf.fast_poly_multiplication(a, b, c, 1 << 3);

	//gf.print_poly(c);


	//d[0] = 1;
	//d[4] = 1;
	//d[6] = 1;
	//d[7] = 1;
	//d[12] = 1;
	//d[13] = 1;
	//d[16] = 1;

	//gf.add_subpoly_with_modular_shift(d, d, 0, 4, 8, 2, 5);


	//gf.SCHONHAGE_DFT(b, c, 3, 3, 3, 6, 0);
	//gf.SCHONHAGE_DFT(c, e, 3, 6, 3, 6, 0);
	//gf.SCHONHAGE_CONVOLUTION(a, b, c, 3, 3, 3, 0);
	//gf.add_subpoly_with_modular_shift(d, d, 0, 3, 6, 3, );

	//gf.fast_poly_multiplication(a, b, c, 4);
	//gf.print_poly(b);
	//gf.print_poly(c);

	//gf.fast_poly_multiplication(a, b, e);
	//gf.print_poly(e);
	//gf.print_poly(e);
	//gf.print_poly(a);

	return 0;
}
//a: 38, 45, 31, 40, 3, 55, 51, 6, 50, 30, 31, 62, 5, 44, 38, 42, 18, 33, 47, 30, 39, 11, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
//b: 56 0 60 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
//3 2
//got: 24 46 25 47 1 19 58 23 49 41 35 31 23 28 59 32 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
//expected: 50 48 25 47 1 19 58 23 49 3 61 31 23 28 59 32 13 44 41 19 51 52 16 46 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
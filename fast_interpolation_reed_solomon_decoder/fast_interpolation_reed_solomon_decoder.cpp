// fast_interpolation_reed_solomon_decoder.cpp : Defines the entry point for the application.
//

#include "fast_interpolation_reed_solomon_decoder.h"
#include "galois_field.h"
#include "fft.h"
#include "berlecamp_massey_decoder.h"
#include <cassert>
#include <bit>
#include <fstream>

InterpolationBasedFastRSDecoder::InterpolationBasedFastRSDecoder(galois_field && gf, unsigned n, unsigned k)
	: _gf(std::move(gf))
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
	//_gf.print_poly(cw);
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
	std::cout << "errors amount: " << t << "\n";
	//assert(t == _t);
	/*if (t < _t) {
		_gf.print_poly(_tmp[0]);
		_gf.print_poly(_tmp[1]);
		_gf.print_poly(_emgcd_tmp[0].first);
		_gf.print_poly(_emgcd_tmp[0].second);
		_gf.print_poly(_emgcd_tmp[1].second);
		_gf.print_poly(_emgcd_tmp[2].second);
	}*/
	std::fill(_tmp[3].begin(), _tmp[3].end(), 0);
	std::fill(_tmp[2].begin(), _tmp[2].end(), 0);
	//std::fill(_tmp[2].begin(), _tmp[2].end(), 0);
	std::fill(_tmp[1].begin(), _tmp[1].end(), 0);

	if (t != 0) {
		std::copy(_tmp[0].begin() + _n - 2 * t + 1, _tmp[0].begin() + _n, _tmp[1].begin());
		std::copy(_tmp[0].begin() + _n - 2 * t, _tmp[0].begin() + _n - t, _tmp[2].begin());
		_gf.SOLVE_TOEPITZ(_tmp[1], _tmp[2], t - 1, _tmp[3]);

		// tmp[3] <=> eta
		std::copy(_tmp[3].begin(), _tmp[3].begin() + t, _tmp[4].begin());
		std::copy(_tmp[0].begin() + _k, _tmp[0].begin() + _k + t, _tmp[5].begin());
		std::reverse(_tmp[4].begin(), _tmp[4].begin() + t + 1); // D(x)
		//std::cout << "before reverse:\n";
		//std::cout << "k: " << _k << " t: " << _t << "\n";
		//_gf.print_poly(_tmp[0]);
		//_gf.print_poly(_tmp[5]);
		std::reverse(_tmp[5].begin(), _tmp[5].begin() + t); // h0 .. h_(t-1)
		_tmp[4][0] = 1;
		_gf.inv_poly(_tmp[4], _tmp[6], _k + t);
		_gf.remainder_of_power(_gf.fast_poly_multiplication(_tmp[5], _tmp[4], _tmp[8]), t);
		_gf.remainder_of_power(_gf.fast_poly_multiplication(_tmp[6], _tmp[8], _tmp[7]), _k + t);
		std::reverse(_tmp[7].begin(), _tmp[7].begin() + t + _k);
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


template<typename DecType>
void test_decoder(DecType & decoder, unsigned n, unsigned k, unsigned t, unsigned iters, std::ofstream& metrics_dump) {
	//encoding::bch_decoder decoder(gf, n, k, 2, n - k + 1, 1);
	//InterpolationBasedFastRSDecoder decoder(gf, n, k);
	//for (size_t t = (n - k ) / 2 - 10; t <= (n - k)/2; ++t) {
		std::cout << "testing t: " << t << "\n";
		//decoder._t = t;
		for (size_t _ = 0; _ < iters; ++_) {
			auto msg = generate_message(n, k);
			std::vector<unsigned> encoded = decoder.encode(msg);
			//gf.DFT(msg, encoded);
			//auto encoded = decoder.encode(msg);
			auto errors = generate_errors(n, t);
			std::vector<unsigned> msg_with_errors(n);
			decoder._gf.add_poly(encoded, errors, msg_with_errors, 0);
			//size_t sz = n * gf._m;
			//linalg::bit_vector vt(sz, false);
			decoder._gf.reset_counters();
			decoder.decode(msg_with_errors);
			//gf.translate_to_bit_vector(msg_with_errors, vt);
			//decoder._decoding.cp(vt);
			//auto res = decoder.decode();
			//gf.translate_bit_vector(res, msg_with_errors);

			std::cout << "additions: " << decoder._gf._additions << " multiplications: " << decoder._gf._multiplications << "\n";
			metrics_dump << n << " " << k << " " << t << " " << decoder._gf._additions << " " << decoder._gf._multiplications << "\n";
			decoder._gf.reset_counters();
			for (size_t i = 0; i < n; ++i) {
				if (encoded[i] != msg_with_errors[i]) {
					std::cout << "decoding error occured with t:"<< t << "\n";
					decoder._gf.print_poly(msg);
					decoder._gf.print_poly(encoded);
					decoder._gf.print_poly(errors);
					decoder._gf.print_poly(msg_with_errors);
					return;
				}
			}
		//}
		//break;
	}
	std::cout << "tests passed\n";


}

// 1: output
// 2: code rate
// 3: error rate
// 4: decoder type: 0 - interpolation decoder, 1 - berlecamp-massey decoder
int main(int argc, char* argv[])
{	
	double k_rate, err_rate;
	int type;
	std::ofstream out("interpolation_decoder_results_0_25rate_max_errs .txt");
	//out << "hello\n";
	if (false) {
		if (argc != 4) {
			std::cout << "expected 3 args: output file, code rate (from 0 to 1), error rate (from 0 to 1)\n";
			return -1;
		}
		try {
			k_rate = std::stod(argv[2]);
			err_rate = std::stod(argv[3]);
			type = std::stoi(argv[4]);
		}
		catch (std::invalid_argument const& e) {
			std::cout << "could not parse numeric args: " << e.what();
			return -1;
		}
	}
	else {
		k_rate = 0.5;
		err_rate = 0.9;
		type = 1;
	}
	if (out.bad()) {
		std::cout << "could not open input file\n";
		return -1;
	}

	out << "n, k, err num, additions, multiplications\n";
	//galois_field gf2(7, 0x83, 7);

	//test_decoder(gf2, 127, 11, 15);

	//galois_field gf(3, 0xb, 3);
	//std::vector<unsigned> field_sizes = { 3, 6, 7, 10 };
	//std::vector<unsigned> field_generators = { 0xb, 0x43, 0x83, 0x2011 };
	//// 2 2
	//for (size_t i = 0; i < field_sizes.size(); ++i) {
	//	galois_field gf(field_sizes[i], field_generators[i], field_sizes[i]);
	//	unsigned n = (1 << field_sizes[i]) - 1;
	//	unsigned k = n * k_rate;
	//	unsigned t = ((n - k) / 2) * err_rate;
	//	if (type == 0) {
	//		InterpolationBasedFastRSDecoder decoder(std::move(gf), n, k);
	//		test_decoder(decoder, n, k, t, 1, out);
	//	}
	//	else if (type == 1) {
	//		encoding::bch_decoder decoder(std::move(gf), n, k, 2, n - k + 1, 1);
	//		std::cout << "here\n";
	//		test_decoder(decoder, n, k, t, 1, out);
	//		//test_decoder(decoder, 511, 256, t, 10, out);
	//	}
	//	else {
	//		std::cout << "unknown decoder type\n";
	//		return -1;
	//	}
	//}
	//InterpolationBasedFastRSDecoder decoder(gf, 511, 256);
	//test_decoder(gf, 7, 1, 10);
	//test_decoder(decoder, 511, 256, 10);

	std::vector<unsigned> a(1024), b(1024), c(80), d(80), e(80), f(80);
	for (size_t i = 1024; i < 2048; ++i) {
		galois_field gf(10, 2011, 10);
		a[1] = 1;
		gf.DFT(a, b);
		std::cout << i << '\n';
		if (b[2] == 4 && b[1] == 2) {
			gf.print_poly(b);
			gf.print_poly(gf._exp_table);
			break;
		}

	}
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

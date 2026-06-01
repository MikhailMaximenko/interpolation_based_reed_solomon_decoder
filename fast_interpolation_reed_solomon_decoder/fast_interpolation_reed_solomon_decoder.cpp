// fast_interpolation_reed_solomon_decoder.cpp : Defines the entry point for the application.
//

#include "fast_interpolation_reed_solomon_decoder.h"
#include "galois_field.h"

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

	std::reverse(_tmp[1].begin(), _tmp[1].begin() + 2 * max_t);
	_tmp[2][2 * max_t] = 1;
	std::cout << "here1\n";
	_gf.EMGCD(_tmp[2], _tmp[1], _emgcd_tmp, 0);
	std::cout << "here2\n";

	size_t fst_deg = _gf.degree(_emgcd_tmp[0].first);
	size_t t = _gf.degree(_emgcd_tmp[2].second);
	std::cout << "errors amount: " << t << "\n";
	std::fill(_tmp[3].begin(), _tmp[3].end(), 0);
	std::fill(_tmp[2].begin(), _tmp[2].end(), 0);
	//std::fill(_tmp[2].begin(), _tmp[2].end(), 0);
	std::fill(_tmp[1].begin(), _tmp[1].end(), 0);

	if (t != 0) {
		std::copy(_tmp[0].begin() + _n - 2 * t + 1, _tmp[0].begin() + _n, _tmp[1].begin());
		std::copy(_tmp[0].begin() + _n - 2 * t, _tmp[0].begin() + _n - t, _tmp[2].begin());
		_gf.SOLVE_TOEPITZ(_tmp[1], _tmp[2], t - 1, _tmp[3]);

		std::copy(_tmp[3].begin(), _tmp[3].begin() + t, _tmp[4].begin());
		std::copy(_tmp[0].begin() + _k, _tmp[0].begin() + _k + t, _tmp[5].begin());
		std::reverse(_tmp[4].begin(), _tmp[4].begin() + t + 1); // D(x)
		std::reverse(_tmp[5].begin(), _tmp[5].begin() + t); // h0 .. h_(t-1)
		_tmp[4][0] = 1;
		std::cout << "operations before inversion:\n";
		std::cout << _gf._additions << " " << _gf._multiplications << "\n";
		_gf.inv_poly(_tmp[4], _tmp[6], _k + t);
		_gf.remainder_of_power(_gf.fast_poly_multiplication(_tmp[5], _tmp[4], _tmp[8]), t);
		_gf.remainder_of_power(_gf.fast_poly_multiplication(_tmp[6], _tmp[8], _tmp[7]), _k + t);
		std::reverse(_tmp[7].begin(), _tmp[7].begin() + t + _k);
		std::copy(_tmp[7].begin(), _tmp[7].begin() + _k, _tmp[0].begin());


		using std::swap;
		swap(cw, _tmp[1]);
		_gf.DFT(_tmp[0], _tmp[2]);
		_gf.add_poly(_tmp[1], _tmp[2], cw, 0);
	}



}

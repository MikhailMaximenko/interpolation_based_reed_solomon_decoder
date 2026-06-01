#include "galois_field.h"
#include "fft.h"
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <vector>
#include <bit>
#include <cassert>

#include <iostream>

galois_field::galois_field(unsigned m, unsigned gen_poly, unsigned poly_size)
	: _m(m)
	, _q(1 << _m)
	, _n(_q - 1)
	, _gen_poly(gen_poly)
	, _poly_size(poly_size)
	, _log_table(_q)
	, _exp_table(_q)
	, _inverse_temporary1(_q)
	, _inverse_temporary2(_q)
	, _inverse_element(_n)
	, _division_tmp(_q)
	, _const2(_q, 0)
	, _a_tmp(_q)
	, _b_tmp(_q)
	, _schonhage_convolution_tmp(5)
	, _schonhage_dft_tmp(5)
	, _schonhage_dft_results_tmp(5)
	, _schonhage_strassen_tmp(5)
	, _multiplication_result_tmp(6 * _n)
	, _gcd_tmp_poly(_q)
	, _caratsuba_tmp(16)
	, _taylor_expansion_tmp(_m + 1)
	, precomputed_basises_delta(_m + 1)
	, precomputed_basises_gamma(_m)
	, precomputed_space_gamma(_m)
	, _gao_mateer_fft_tmp(_m + 1)
{
	init();
}

void galois_field::init() {
	_const2[0] = 2;
	unsigned x = 1;
	_exp_table[_n] = 1;
	//_log_table[0] = 0;
	for (int i = 0; i < _n; i++) {
		_exp_table[i] = x;
		_log_table[x] = i;

		// Multiply by generator (x)
		x <<= 1;
		x = x ^ ((x & (1 << _poly_size)) ? _gen_poly : 0);
	}

	// code for cantor basis
	// evaluate s
	// recurrent way of evaluation may be easier
	//_s.resize(_m + 1);
	//_s.back().resize(_q);
	//_s[0].resize(_q);
	//for (size_t i = 1; i < _q; ++i) {
	//	_s[0][i] = i;
	//}
	//for (size_t level = 1; level < _m; ++level) {
	//	_s[level].resize(_q);
	//	for (size_t i = 0; i < _q; ++i) {
	//		_s[level][i] = 1;
	//		for (size_t k = 0; k < (1 << level); ++k) {
	//			_s[level][i] = multiply(_s[level][i], add(k, i));
	//		}
	//	}
	//}
	////std::cout << "s:" << _s.size() << "\n";
	//for (size_t level = 0; level < _m; ++level) {
	//	//std::cout << level << ": ";
	//	unsigned inv = inverse(_s[level][1 << level]);
	//	for (size_t i = 0; i < _q; ++i) {
	//		//_s[level][i] = multiply(_s[level][i], inv);
	//		//std::cout << _s[level][i] << " ";
	//	}
	//	//std::cout << "\n";
	//}

	////std::cout << 3 << ": ";
	//for (size_t i = 0; i < _q; ++i) {
	//	//std::cout << _s[3][i] << " ";
	//}
	////std::cout << "\n";
	size_t tmp_sizes = _q;
	_emgcd_tmp_polynomials.resize(20);
	for (auto& tmp : _emgcd_tmp_polynomials) {
		for (auto& v : tmp) {
			v.resize(tmp_sizes);
		}
	}

	_emgcd_tmp_result.resize(20);
	for (auto& tmp : _emgcd_tmp_result) {
		for (auto& v : tmp) {
			v.first.resize(tmp_sizes);
			v.second.resize(tmp_sizes);
		}
	}
	_emgcd_tmp_result2.resize(20);
	for (auto& tmp : _emgcd_tmp_result2) {
		for (auto& v : tmp) {
			v.first.resize(tmp_sizes);
			v.second.resize(tmp_sizes);
		}
	}

	_dft_tmp.resize(_m + 1);
	_idft_tmp.resize(_m + 1);
	for (auto& a : _dft_tmp) {
		a.resize(tmp_sizes);
		for (auto& v : a) {
			v.resize(tmp_sizes);
		}
	}
	for (auto& a : _idft_tmp) {
		a.resize(tmp_sizes);
		for (auto& v : a) {
			v.resize(tmp_sizes);
		}
	}

	for (auto& v : _ad_tmp_polinomyals) {
		v.resize(tmp_sizes);
	}
	for (auto& tmp : _ad_tmp_emgcd_results) {
		for (auto& v : tmp) {
			v.first.resize(tmp_sizes);
			v.second.resize(tmp_sizes);
		}
	}

	for (auto& v : _gcd_tmp) {
		v.first.resize(3 * _q);
		v.second.resize(3 * _q);
	}
	for (auto& tmp : _solve_toeplitz_tmp) {
		tmp.resize(tmp_sizes);
	}

	std::cout << "initting mult tmps\n";
	size_t size = 9 * _q;

	for (size_t i = 0; i < _schonhage_dft_tmp.size(); ++i) {
		_schonhage_dft_tmp[i].resize(size);
		for (auto& v : _schonhage_dft_tmp[i]) {
			v.resize(size);
		}
		for (auto& v : _schonhage_convolution_tmp[i]) {
			v.resize(size);
		}
		_schonhage_dft_results_tmp[i].resize(10);
		for (auto& v : _schonhage_dft_results_tmp[i]) {
			v.resize(size);
		}
		for (auto& v : _schonhage_strassen_tmp[i]) {
			v.resize(size);
		}
		size /= 3;
		++size;
	}
	size = 1 << (sizeof(unsigned) * 8 - std::countl_zero(_n) + 4);
	for (size_t i = 0; i < _caratsuba_tmp.size(); ++i) {
		for (auto& v : _caratsuba_tmp[i]) {
			v.resize(size);
		}
		size >>= 1;
	}

	std::cout << "initting gao-mateer subspaces\n";
	unsigned max_m = _m;
	for (size_t i = max_m; i >= 1; --i) {
		for (auto& v : _taylor_expansion_tmp[i]) {
			v.resize(1 << i);
		}
	}
	for (size_t i = max_m; i >= 1; --i) {
		for (auto& v : _gao_mateer_fft_tmp[i]) {
			v.resize(1 << i);
		}
	}
	// init gao-mateer
	// init the biggest space
	precomputed_basises_delta[max_m].resize(max_m);
	for (size_t i = 0; i < max_m; ++i) {
		precomputed_basises_delta[max_m][i] = (1 << i);
	}
	for (size_t m = max_m - 1; m >= 1; --m) {
		precomputed_basises_delta[m].resize(m);
		precomputed_basises_gamma[m].resize(m);
		precomputed_space_gamma[m].resize(1 << m);
		for (size_t i = 0; i < m; ++i) {
			precomputed_basises_gamma[m][i] = multiply(precomputed_basises_delta[m + 1][i], inverse(precomputed_basises_delta[m + 1][m]));
			precomputed_basises_delta[m][i] = add(multiply(precomputed_basises_gamma[m][i], precomputed_basises_gamma[m][i]), precomputed_basises_gamma[m][i]);

			for (size_t j = 0; j < (1 << m); ++j) {
				if ((j & (1 << i)) != 0) {
					precomputed_space_gamma[m][j] = add(precomputed_space_gamma[m][j], precomputed_basises_gamma[m][i]);
				}

			}
		}
	}

	std::cout << "initted\n";
}

// legacy code for arbitrary fields

//std::vector<unsigned> galois_field::inverse_add_init(std::vector<unsigned> const& v) const {
//	std::vector<unsigned> res;
//	for (unsigned i = 0; i < v.size(); ++i) {
//		res.emplace_back((_p - v[i]) % _p);
//	}
//	return  res;
//}

//std::vector<unsigned> galois_field::add_init(std::vector<unsigned> const& a, std::vector<unsigned> const& b) const {
//	std::vector<unsigned> res;
//	unsigned i = 0;
//	while (i < a.size() && i < b.size()) {
//		res.emplace_back((a[i] + b[i]) % _p);
//		++i;
//	}
//	while (i < a.size()) {
//		res.emplace_back(a[i]);
//		++i;
//	}
//	while (i < b.size()) {
//		res.emplace_back(b[i]);
//		++i;
//	}
//	return res;
//}

//std::vector<unsigned> galois_field::multiply_by_const_init(std::vector<unsigned> const& a, unsigned b) const {
//	std::vector<unsigned> res;
//	for (unsigned i : a) {
//		res.push_back(i * b % _p);
//	}
//	return res;
//}

unsigned galois_field::degree(std::vector<unsigned> const& a) {
	return degree(a, 0, a.size());
}

unsigned galois_field::degree(std::vector<unsigned> const& a, unsigned start, unsigned end) {
	if (end <= start || a.size() == 0) { return 0; }
	for (unsigned deg = std::min((size_t)end, a.size()) - 1; deg > start; --deg) {
		if (a[deg] != 0) {
			return deg - start;
		}
	}
	return 0;
}

unsigned galois_field::poly_to_num(std::vector<unsigned> const& a) const {
	unsigned res = 0;
	for (ptrdiff_t i = a.size() - 1; i >= 0; --i) {
		res *= 2;
		res += a[i];
	}
	return res;
}

std::vector<unsigned>& galois_field::shift_poly(std::vector<unsigned>& a) const {
	for (ptrdiff_t i = a.size() - 1; i >= 1; --i) {
		a[i] = a[i - 1];
	}
	a[0] = 0;
	return a;
}

unsigned galois_field::add(unsigned a, unsigned b) const { 
	++_additions;
	return a ^ b;
}

unsigned galois_field::multiply(unsigned a, unsigned b) const { 
	++_multiplications;
	if (a == 0 || b == 0) {
		return 0;
	}
	unsigned log_sum = (_log_table[a] + _log_table[b]) % _n;
	return _exp_table[log_sum]; 
}

unsigned galois_field::multiplyConst(unsigned a, int x) const {
	++_multiplications;
	if (x < 0) {
		return 0;
	}
	return a ? _exp_table[(_log_table[a] + x) % _n] : 0;
}
unsigned galois_field::inverse(unsigned a) const {
	assert(a != 0);
	if (a == 0) throw std::runtime_error("division by zero");
    return _exp_table[_n - _log_table[a]];
}
unsigned galois_field::divide(unsigned a, unsigned b) const { 
	return multiply(a, inverse(b)); 
}

// deprecated
std::vector<unsigned>& galois_field::fast_poly_multiplication(std::vector<unsigned>& a, std::vector<unsigned>& b, std::vector<unsigned>& dst) {

	auto dg = std::max(degree(a), degree(b)) + 1;

	size_t answer_dg = degree(a) + degree(b);
	size_t start_add = _additions;
	size_t start_mult = _multiplications;
	fast_poly_multiplication(a, b, dst, dg);

	std::fill(dst.begin() + std::min(answer_dg + 1, dst.size()), dst.end(), 0);
	
	return dst;
}

std::vector<unsigned>& galois_field::fast_poly_multiplication(std::vector<unsigned>& a, std::vector<unsigned>& b, std::vector<unsigned>& dst, unsigned length) {
	++_poly_multiplications;
	if (length <= 256) {
		size_t len = 1 << (sizeof(unsigned) * 8 - std::countl_zero(length - 1));
		caratsuba_multiplication(a, b, _multiplication_result_tmp, len, 0);

	}
	else {
		unsigned len = length - 1;
		unsigned n = 1;
		while (len >= 1) {
			n *= 3;
			len /= 3;
		}


		SCHONHAGE_STRASSEN_FFT(a, b, _multiplication_result_tmp, n, 0);
	}
	/*if ((length <= 512 && (length > _n / 2 || length < _n / 4)) || length >= 64) {
	}
	else if (length <= _q / 2 && length >= _q / 4) {
		for (size_t i = 0; i < 3; ++i) {
			std::fill(_caratsuba_tmp[0][i].begin(), _caratsuba_tmp[0][i].begin() + _n, 0);
		}
		std::fill(_multiplication_result_tmp.begin(), _multiplication_result_tmp.end(), 0);
		exDFT(a, _caratsuba_tmp[0][0]);
		exDFT(b, _caratsuba_tmp[0][1]);
		for (size_t i = 0; i < _q; ++i) {
			_caratsuba_tmp[0][2][i] = multiply(_caratsuba_tmp[0][0][i], _caratsuba_tmp[0][1][i]);
		}

		exIDFT(_caratsuba_tmp[0][2], _multiplication_result_tmp);
		for (size_t i = 0; i < 3; ++i) {
			std::fill(_caratsuba_tmp[0][i].begin(), _caratsuba_tmp[0][i].begin() + _n, 0);
		}
		
	
	} else {
		unsigned len = length - 1;
		unsigned n = 1;
		while (len >= 1) {
			n *= 3;
			len /= 3;
		}


		SCHONHAGE_STRASSEN_FFT(a, b, _multiplication_result_tmp, n, 0);
	}*/
	std::copy(_multiplication_result_tmp.begin(), _multiplication_result_tmp.begin() + std::min((size_t)2ull * length, dst.size()), dst.begin());
	return dst;
}

std::vector<unsigned>& galois_field::fast_poly_division(std::vector<unsigned>& a, std::vector<unsigned>& b, std::vector<unsigned>& quotient, std::vector<unsigned>& remainder) {
	++_poly_divisions;
	auto deg_a = degree(a), deg_b = degree(b);
	std::fill(remainder.begin(), remainder.end(), 0);
	std::fill(quotient.begin(), quotient.end(), 0);
	if (deg_a < deg_b) {
		std::copy(a.begin(), a.end(), remainder.begin());

		return a;
	}
	auto m = deg_a - deg_b;

	rev_poly(a, quotient, deg_a);
	rev_poly(b, _division_tmp, deg_b);
	inv_poly(_division_tmp, remainder, m + 1);
	fast_poly_multiplication(quotient, remainder, _division_tmp);
	remainder_of_power(_division_tmp, m + 1);
	rev_poly(_division_tmp, quotient, m);
	fast_poly_multiplication(quotient, b, _division_tmp);
	add_poly(a, _division_tmp, remainder, 0);
	return quotient;
}

void galois_field::reset_counters() {
	std::cout << "resetting with: \n\tadditions: " << _additions << "\n\tmultiplications: "
		<< _multiplications << "\n\tpoly multiplications: " << _poly_multiplications << "\n\tpoly divisions: " << _poly_divisions << "\n";
	_additions = 0;
	_multiplications = 0;
	_poly_divisions = 0;
	_poly_multiplications = 0;
}

std::array<std::pair<std::vector<unsigned>, std::vector<unsigned>>, 3>& 
		galois_field::EMGCD(std::vector<unsigned> const& u0, std::vector<unsigned> const& u1, std::array<std::pair<std::vector<unsigned>, std::vector<unsigned>>, 3>& dst, unsigned tmp_num) {
	unsigned n = degree(u0);
	auto deg_u1 = degree(u1);
	//std::cout << "call " << tmp_num << "\n";
	//print_poly(u0);
	//print_poly(u1);
	unsigned m = n / 2 + n % 2;

	if ((deg_u1 == 0 && u1[0] == 0) || deg_u1 < m) {
		dst[0].first = u0;
		dst[0].second = u1;
		std::fill(dst[1].first.begin(), dst[1].first.end(), 0);
		std::fill(dst[2].first.begin(), dst[2].first.end(), 0);
		std::fill(dst[1].second.begin(), dst[1].second.end(), 0);
		std::fill(dst[2].second.begin(), dst[2].second.end(), 0);

		dst[1].first[0] = 1;
		dst[2].second[0] = 1;
		//std::cout << "ret here\n";
		return dst;
	}
	auto& locals = _emgcd_tmp_polynomials[tmp_num];
	for (auto& v : locals) {
		std::fill(v.begin(), v.end(), 0);
	}
	std::array<std::pair<std::vector<unsigned>, std::vector<unsigned>>, 3>& result = _emgcd_tmp_result[tmp_num];
	// tmp layout:
	// b0   b1   c0   c1   d    e    q    f    
	// 0    1    2    3    4    5    6    7  
	
	split_poly(u0,locals[0], locals[2], m);
	split_poly(u1, locals[1], locals[3], m);
	
	// u1w1v1 is stored in _emgcd_tmp_result
	EMGCD(locals[0], locals[1], _emgcd_tmp_result[tmp_num], tmp_num + 1);

	using std::swap;
	
	// b0, b1 wont be used further
	fast_poly_multiplication(result[1].first, locals[2], locals[4]);
	fast_poly_multiplication(result[2].first, locals[3], locals[1]);
	
	add_poly(locals[4], locals[1], locals[0], 0);
	add_poly(locals[0], result[0].first, locals[4], m); // d
	
	// in some locals may be litter
	fast_poly_multiplication(result[1].second, locals[2], locals[5]);
	fast_poly_multiplication(result[2].second, locals[3], locals[1]);
	add_poly(locals[5], locals[1], locals[0], 0);
	add_poly(locals[0], result[0].second, locals[5], m); // e

	unsigned deg_e = degree(locals[5]);
	if (deg_e < m) {
		dst[0].first = locals[4];
		dst[0].second = locals[5];
		dst[1].first = result[1].first;
		dst[1].second = result[1].second;
		dst[2].first = result[2].first;
		dst[2].second = result[2].second;
		return dst;
	}

	// c0 and c1 also can be reused -> g0,g1,h0,h1 will be stored there (in b0,b1,c0,c1)

	// tmp layout:
	// g0   g1   h0   h1   d    e    q    f    
	// 0    1    2    3    4    5    6    7  
	fast_poly_division(locals[4], locals[5], locals[6], locals[7]);
	unsigned k = 2 * m - deg_e;
	split_poly(locals[5], locals[0], locals[2], k);
	split_poly(locals[7], locals[1], locals[3], k);
	EMGCD(locals[0], locals[1], _emgcd_tmp_result2[tmp_num], tmp_num+1); // need to be careful with such storage

	std::array<std::pair<std::vector<unsigned>, std::vector<unsigned>>, 3>& result2 = _emgcd_tmp_result2[tmp_num];
	//g0, g1 wont be used further
	
	// fst part of the answer
	fast_poly_multiplication(result2[1].first, locals[2], locals[0]);
	fast_poly_multiplication(result2[2].first, locals[3], locals[1]);
	add_poly(locals[0], locals[1], locals[4], 0);
	swap(locals[0], locals[4]);

	add_poly(locals[0], result2[0].first, locals[4], k); // d

	// in some locals may be litter
	fast_poly_multiplication(result2[1].second, locals[2], locals[0]);
	fast_poly_multiplication(result2[2].second, locals[3], locals[1]);
	add_poly(locals[0], locals[1], locals[5], 0);
	swap(locals[0], locals[5]);

	add_poly(locals[0], result2[0].second, locals[5], k); // e
	// locals [0..3] can be reused

	dst[0].first = locals[4];
	dst[0].second = locals[5];

	// multiply matrices (may be optimized)
	fast_poly_multiplication(result2[2].first, locals[6], locals[0]);
	add_poly(result2[1].first, locals[0], locals[1], 0);

	fast_poly_multiplication(result2[2].second, locals[6], locals[0]);
	add_poly(result2[1].second, locals[0], locals[2], 0);

	fast_poly_multiplication(result2[2].first, result[1].first, locals[0]);
	fast_poly_multiplication(locals[1], result[1].second, locals[3]);
	add_poly(locals[0], locals[3], dst[1].first, 0);

	fast_poly_multiplication(result2[2].second, result[1].first, locals[0]);
	fast_poly_multiplication(locals[2], result[1].second, locals[3]);
	add_poly(locals[0], locals[3], dst[1].second, 0);

	fast_poly_multiplication(result2[2].first, result[2].first, locals[0]);
	fast_poly_multiplication(locals[1], result[2].second, locals[3]);
	add_poly(locals[0], locals[3], dst[2].first, 0);

	fast_poly_multiplication(result2[2].second, result[2].first, locals[0]);
	fast_poly_multiplication(locals[2], result[2].second, locals[3]);
	add_poly(locals[0], locals[3], dst[2].second, 0);
	// end multiply 


	return dst;
}

std::vector<unsigned>& galois_field::GCD(std::vector<unsigned> const& a, std::vector<unsigned> const& b, std::vector<unsigned>& dst) {

	size_t a_deg = degree(a);
	size_t b_deg = degree(b);
	if (a_deg >= b_deg) {
		std::copy(a.begin(), a.begin() + degree(a) + 1, dst.begin());
		std::copy(b.begin(), b.begin() + degree(b) + 1, _gcd_tmp_poly.begin());
	}
	else {
		std::copy(a.begin(), a.begin() + degree(a) + 1, _gcd_tmp_poly.begin());
		std::copy(b.begin(), b.begin() + degree(b) + 1, dst.begin());
		std::swap(a_deg, b_deg);
	}
	do {
		if (b_deg * 2 <= a_deg || b_deg == a_deg) {
			fast_poly_division(dst, _gcd_tmp_poly, _gcd_tmp[0].first, _gcd_tmp[0].second);
			std::copy(_gcd_tmp_poly.begin(), _gcd_tmp_poly.begin() + a_deg + 1, dst.begin());
			std::copy(_gcd_tmp[0].second.begin(), _gcd_tmp[0].second.begin() + b_deg + 1, _gcd_tmp_poly.begin());
		}
		else {
			EMGCD(dst, _gcd_tmp_poly, _gcd_tmp, 0);
			std::copy(_gcd_tmp[0].first.begin(), _gcd_tmp[0].first.begin() + a_deg + 1, dst.begin());
			std::copy(_gcd_tmp[0].second.begin(), _gcd_tmp[0].second.begin() + b_deg + 1, _gcd_tmp_poly.begin());
		}
		a_deg = degree(dst);
		b_deg = degree(_gcd_tmp_poly);
	} while (b_deg != 0 || _gcd_tmp_poly[0] != 0);

	return dst;
}

void galois_field::AD(std::vector<unsigned>& a, unsigned n, std::vector<unsigned>& dst1, std::vector<unsigned>& dst2) {
	for (auto& v : _ad_tmp_polinomyals) {
		std::fill(v.begin(), v.end(), 0);
	}
	_ad_tmp_polinomyals[0][2 * n + 1] = 1;
	_ad_tmp_polinomyals[1] = a;
	EMGCD(_ad_tmp_polinomyals[0], _ad_tmp_polinomyals[1], _ad_tmp_emgcd_results[0], 0);
	if (degree(_ad_tmp_emgcd_results[0][0].second) < n) {
		throw std::runtime_error("singular");
	}
	if (_ad_tmp_emgcd_results[0][2].second[0] != 0) {
		rev_poly(_ad_tmp_polinomyals[1], _ad_tmp_polinomyals[2], 2 * n);
		EMGCD(_ad_tmp_polinomyals[0], _ad_tmp_polinomyals[2], _ad_tmp_emgcd_results[1], 0);
		dst1 = multipy_poly_by_const(_ad_tmp_emgcd_results[0][2].second, inverse(_ad_tmp_emgcd_results[0][0].second[degree(_ad_tmp_emgcd_results[0][0].second)]));
		dst2 = multipy_poly_by_const(_ad_tmp_emgcd_results[1][2].second, inverse(_ad_tmp_emgcd_results[1][0].second[degree(_ad_tmp_emgcd_results[1][0].second)]));
		return;
	}
	_ad_tmp_polinomyals[0][2 * n + 1] = 0;
	_ad_tmp_polinomyals[0][2 * n + 3] = 1;
	shift_poly(_ad_tmp_polinomyals[1]);
	_ad_tmp_polinomyals[1][2 * n + 2] = 1;
	_ad_tmp_polinomyals[1][0] = 1;
	EMGCD(_ad_tmp_polinomyals[0], _ad_tmp_polinomyals[1], _ad_tmp_emgcd_results[0], 0);
	if (degree(_ad_tmp_emgcd_results[0][0].second) == n + 1) {
		rev_poly(_ad_tmp_polinomyals[1], _ad_tmp_polinomyals[2], 2 * n + 2);
		EMGCD(_ad_tmp_polinomyals[0], _ad_tmp_polinomyals[2], _ad_tmp_emgcd_results[1], 0);
		dst1 = multipy_poly_by_const(_ad_tmp_emgcd_results[0][2].second, inverse(_ad_tmp_emgcd_results[0][0].second[degree(_ad_tmp_emgcd_results[0][0].second)]));
		dst2 = multipy_poly_by_const(_ad_tmp_emgcd_results[1][2].second, inverse(_ad_tmp_emgcd_results[1][0].second[degree(_ad_tmp_emgcd_results[1][0].second)]));
		return;
	}
	_ad_tmp_polinomyals[1][2 * n + 2] = _inverse_element;
	EMGCD(_ad_tmp_polinomyals[0], _ad_tmp_polinomyals[1], _ad_tmp_emgcd_results[0], 0);
	rev_poly(_ad_tmp_polinomyals[1], _ad_tmp_polinomyals[2], 2 * n + 2);
	EMGCD(_ad_tmp_polinomyals[0], _ad_tmp_polinomyals[2], _ad_tmp_emgcd_results[1], 0);
	dst1 = multipy_poly_by_const(_ad_tmp_emgcd_results[0][2].second, inverse(_ad_tmp_emgcd_results[0][0].second[degree(_ad_tmp_emgcd_results[0][0].second)]));
	dst2 = multipy_poly_by_const(_ad_tmp_emgcd_results[1][2].second, inverse(_ad_tmp_emgcd_results[1][0].second[degree(_ad_tmp_emgcd_results[1][0].second)]));
}

std::vector<unsigned>& galois_field::SOLVE_TOEPITZ(std::vector<unsigned>& a, std::vector<unsigned>& b, unsigned n, std::vector<unsigned>& dst) {
	AD(a, n, _ad_x, _ad_y);
	for (auto& v : _solve_toeplitz_tmp) {
		std::fill(v.begin(), v.end(), 0);
	}
	assert(_ad_x[0] != 0);
	rev_poly(b, _solve_toeplitz_tmp[0], n); // br
	rev_poly(_ad_x, _solve_toeplitz_tmp[1], n + 1); // RevPoly(xy[1], n + 1)

	remainder_of_power(fast_poly_multiplication(_solve_toeplitz_tmp[1], _solve_toeplitz_tmp[0], _solve_toeplitz_tmp[2]), n + 1); // xn
	remainder_of_power(fast_poly_multiplication(_ad_y, _solve_toeplitz_tmp[0], _solve_toeplitz_tmp[3]), n + 1); // y0

	rev_poly(_solve_toeplitz_tmp[3], _solve_toeplitz_tmp[4], n); // RevPoly(y0, n)
	fast_poly_multiplication(_ad_x, _solve_toeplitz_tmp[4], _solve_toeplitz_tmp[1]); // xy[1]*RevPoly(y0, n)

	rev_poly(_ad_y, _solve_toeplitz_tmp[3], n + 1); // RevPoly(xy[2], n + 1)
	rev_poly(_solve_toeplitz_tmp[2], _solve_toeplitz_tmp[0], n); // RevPoly(xn, n)
	fast_poly_multiplication(_solve_toeplitz_tmp[3], _solve_toeplitz_tmp[0], _solve_toeplitz_tmp[2]); // RevPoly(xy[2], n + 1)*RevPoly(xn, n)
	return multipy_poly_by_const(remainder_of_power(add_poly(_solve_toeplitz_tmp[1], _solve_toeplitz_tmp[2], dst, 0), n + 1), inverse(_ad_x[0]));
}


std::vector<unsigned>& galois_field::DFT(std::vector<unsigned>& src, std::vector<unsigned>& dst) {
	// call fft
	call_fft(*this, src, dst, _n);
	return dst;
}

std::vector<unsigned>& galois_field::IDFT(std::vector<unsigned>& src, std::vector<unsigned>& dst) {
	DFT(src, dst);
	std::reverse(dst.begin() + 1, dst.begin() + _n);
	return dst;
}

std::vector<unsigned>& galois_field::exDFT(std::vector<unsigned>& src, std::vector<unsigned>& dst) {
	std::fill(_dft_tmp[0][0].begin(), _dft_tmp[0][0].end(), 0);
	std::fill(dst.begin(), dst.end(), 0);
	std::copy(src.begin() + 1, src.begin() + _q, _dft_tmp[0][0].begin());
	DFT(_dft_tmp[0][0], dst);
	std::memmove(dst.data() + 1, dst.data(), _n * sizeof(unsigned));
	dst[0] = src[0];
	for (size_t i = 1; i <= _n; ++i) {
		dst[i] = add(multiplyConst(dst[i], i - 1), dst[0]);
	}
	return dst;
}

std::vector<unsigned>& galois_field::exIDFT(std::vector<unsigned>& src, std::vector<unsigned>& dst) {
	std::fill(_dft_tmp[0][0].begin(), _dft_tmp[0][0].end(), 0);
	std::fill(dst.begin(), dst.end(), 0);
	std::copy(src.begin() + 1, src.begin() + _q, _dft_tmp[0][0].begin());
	for (size_t i = 0; i < _n; ++i) {
		_dft_tmp[0][0][i] = multiply(add(_dft_tmp[0][0][i], src[0]), _exp_table[_n - i]);
	}
	IDFT(_dft_tmp[0][0], dst);
	std::memmove(dst.data() + 1, dst.data(), _n * sizeof(unsigned));
	dst[0] = src[0];
	return dst;
}


std::vector<unsigned>& galois_field::DFT(std::vector<unsigned>& src, std::vector<unsigned>& dst, unsigned size, unsigned bias) {
	DFTimpl(src, dst, size, 0, 0);
	return dst;
}

std::vector<unsigned>& galois_field::IDFT(std::vector<unsigned>& src, std::vector<unsigned>& dst, unsigned size, unsigned bias) {
	IDFTimpl(src, dst, size, 0, 0);
	return dst;
}

std::vector<unsigned>& galois_field::add_subpoly(std::vector<unsigned>& a, std::vector<unsigned>& b, unsigned m, unsigned start, unsigned end) {
	for (size_t i = start; i < end; ++i) {
		a[i - start + m] = add(a[i - start + m], b[i]);
	}
	return a;
}

// performs addition part of poly a beginning from m and subpoly of b [start, end) shifted by shift in the ring GF(_q)[x]/x^2mod+x^mod+1
std::vector<unsigned>& galois_field::add_subpoly_with_modular_shift(std::vector<unsigned>& a, std::vector<unsigned>& b, unsigned m, unsigned start, unsigned end, unsigned mod, unsigned shift) {
	if (start >= end) {
		return a;
	}
	shift = shift % (3 * mod);
	for (size_t i = start; i < end; ++i) {
		if (shift < 2 * mod) {
			a[m + shift] = add(a[m + shift], b[i]);
		}
		else {
			a[m + shift - 2 * mod] = add(a[m + shift - 2 * mod], b[i]);
			a[m + shift - mod] = add(a[m + shift - mod], b[i]);
		}
		++shift;
		if (shift == 3 * mod) {
			shift = 0;
		}
	}
	return a;
}

std::vector<unsigned>& galois_field::SCHONHAGE_DFT(std::vector<unsigned>& src, std::vector<unsigned>& dst, unsigned length, unsigned root, unsigned m, unsigned block_size, unsigned lvl) { 
	if (length == 1) {
		std::copy(src.begin(), src.begin() + block_size, dst.begin());
		return dst;
	}
	root %= 3 * m;
	for (auto& v : _schonhage_dft_tmp[lvl]) {
		std::fill(v.begin(), v.end(), 0);
	}
	for (auto& v : _schonhage_dft_results_tmp[lvl]) {
		std::fill(v.begin(), v.end(), 0);
	}
	size_t ln = length / 3;
	// init b and c
	for (size_t i = 0; i < ln; ++i) {
		std::copy(src.begin() + (3 * i) * block_size, src.begin() + (3 * i + 1) * block_size, _schonhage_dft_tmp[lvl][0].begin() + i * block_size);
	}
	SCHONHAGE_DFT(_schonhage_dft_tmp[lvl][0], _schonhage_dft_results_tmp[lvl][0], ln, root * 3, m, block_size, lvl + 1);
	for (size_t i = 0; i < ln; ++i) {
		std::copy(src.begin() + (3 * i + 1) * block_size, src.begin() + (3 * i + 2) * block_size, _schonhage_dft_tmp[lvl][1].begin() + i * block_size);
	}
	SCHONHAGE_DFT(_schonhage_dft_tmp[lvl][1], _schonhage_dft_results_tmp[lvl][1], ln, root * 3, m, block_size, lvl + 1);
	for (size_t i = 0; i < ln; ++i) {
		std::copy(src.begin() + (3 * i + 2) * block_size, src.begin() + (3 * i + 3) * block_size, _schonhage_dft_tmp[lvl][2].begin() + i * block_size);
	}
	SCHONHAGE_DFT(_schonhage_dft_tmp[lvl][2], _schonhage_dft_results_tmp[lvl][2], ln, root * 3, m, block_size, lvl + 1);
	
	unsigned cur = 0;
	for (size_t i = 0; i < ln; ++i) {
		// a[i] =  b[i] + cur*c[i] + cur*cur*d[i]
		// a[i + n/3] =  b[i] + root^(n/3)cur*c[i] + root^(2n/3)cur*cur*d[i]
		// a[i + 2n/3] =  b[i] + root^(2n/3)cur*c[i] + root^(n/3)cur*cur*d[i]

		// x^3m === x^2m + x^m === 1 (mod x^2m + x^m + 1)
		// thus cur := cur % 3m;
		size_t cc = cur;
		size_t ccc = cur * 2 % (3 * m);
		add_subpoly(dst, _schonhage_dft_results_tmp[lvl][0], i * block_size, i * block_size, (i + 1) * block_size);
		add_subpoly_with_modular_shift(dst, _schonhage_dft_results_tmp[lvl][1], i * block_size, i * block_size, (i + 1) * block_size, m, cc);
		add_subpoly_with_modular_shift(dst, _schonhage_dft_results_tmp[lvl][2], i * block_size, i * block_size, (i + 1) * block_size, m, ccc);

		cc += ln * root;
		ccc += (2 * ln) * root;

		add_subpoly(dst, _schonhage_dft_results_tmp[lvl][0], (i + ln) * block_size, i * block_size, (i + 1) * block_size);
		add_subpoly_with_modular_shift(dst, _schonhage_dft_results_tmp[lvl][1], (i + ln) * block_size, i * block_size, (i + 1) * block_size, m, cc);
		add_subpoly_with_modular_shift(dst, _schonhage_dft_results_tmp[lvl][2], (i + ln) * block_size, i * block_size, (i + 1) * block_size, m, ccc);

		cc += ln * root;
		ccc -= ln * root;

		add_subpoly(dst, _schonhage_dft_results_tmp[lvl][0], (i + (2 * ln)) * block_size, i * block_size, (i + 1) * block_size);
		add_subpoly_with_modular_shift(dst, _schonhage_dft_results_tmp[lvl][1], (i + (2 * ln)) * block_size, i * block_size, (i + 1) * block_size, m, cc);
		add_subpoly_with_modular_shift(dst, _schonhage_dft_results_tmp[lvl][2], (i + (2 * ln)) * block_size, i * block_size, (i + 1) * block_size, m, ccc);

		cur += root;
		cur = cur % (3 * m);
	}


	return dst;
}

std::vector<unsigned>& galois_field::SCHONHAGE_CONVOLUTION(std::vector<unsigned>& a, std::vector<unsigned>& b, std::vector<unsigned>& dst, unsigned length, unsigned root, unsigned m, unsigned level) {
	auto& tmp = _schonhage_convolution_tmp[level];
	for (auto& v : tmp) {
		std::fill(v.begin(), v.end(), 0);
	}
	root = root % (3 * m);
	SCHONHAGE_DFT(a, tmp[0], length, root, m, 2 * m, 0);
	SCHONHAGE_DFT(b, tmp[1], length, root, m, 2 * m, 0);
	
	auto block_size = 2 * m;
	for (size_t i = 0; i < length; ++i) {
		std::copy(tmp[0].begin() + i * block_size, tmp[0].begin() + (i + 1) * block_size, tmp[2].begin());
		std::copy(tmp[1].begin() + i * block_size, tmp[1].begin() + (i + 1) * block_size, tmp[3].begin());
		std::fill(tmp[4].begin(), tmp[4].begin() + block_size, 0);
		SCHONHAGE_STRASSEN_FFT(tmp[2], tmp[3], tmp[4], m, level + 1);
		std::copy(tmp[4].begin(), tmp[4].begin() + block_size, tmp[5].begin() + i * block_size);
		add_subpoly_with_modular_shift(tmp[5], tmp[4], 0, block_size, 2 * block_size, m, block_size);
	}
	SCHONHAGE_DFT(tmp[5], dst, length, 3 * m - root, m, 2 * m, 0);
	return dst;
}

std::vector<unsigned>& galois_field::caratsuba_multiplication(std::vector<unsigned>& a, std::vector<unsigned>& b, std::vector<unsigned>& dst, unsigned length, unsigned lvl) {
	auto& tmp = _caratsuba_tmp[lvl];
	std::fill(dst.begin(), dst.end(), 0);
	if (length <= 16) {
		for (size_t i = 0; i < length; ++i) {
			for (size_t j = 0; j < length; ++j) {
				dst[i + j] = add(dst[i + j], multiply(a[i], b[j]));
			}
		}
		return dst;
	}
	//else if (length <= (_q >> 1) && length >= _q / 8) {
	//	for (size_t i = 0; i < 4; ++i) {
	//		std::fill(tmp[i].begin(), tmp[i].begin() + _n, 0);
	//	}
	//	std::fill(_multiplication_result_tmp.begin(), _multiplication_result_tmp.end(), 0);
	//	exDFT(a, tmp[0]);
	//	exDFT(b, tmp[1]);
	//	for (size_t i = 0; i < _q; ++i) {
	//		tmp[2][i] = multiply(tmp[0][i], tmp[1][i]);
	//	}

	//	exIDFT(tmp[2], tmp[3]);
	//	/*for (size_t i = 0; i < 3; ++i) {
	//		std::fill(_caratsuba_tmp[0][i].begin(), _caratsuba_tmp[0][i].begin() + _n, 0);
	//	}*/
	//	std::copy(tmp[3].begin(), tmp[3].begin() + std::min((size_t)length << 1, dst.size()), dst.begin());
	//	return dst;
	//}
	for (auto& v : tmp) {
		std::fill(v.begin(), v.begin() + (length << 1), 0);
	}
	unsigned k = length >> 1;
	std::copy_n(a.begin(), k, tmp[0].begin());
	std::copy_n(a.begin() + k, k, tmp[1].begin());
	std::copy_n(b.begin(), k, tmp[2].begin());
	std::copy_n(b.begin() + k, k, tmp[3].begin());
	caratsuba_multiplication(tmp[0], tmp[2], tmp[4], k, lvl + 1); // F0*G0
	caratsuba_multiplication(tmp[1], tmp[3], tmp[5], k, lvl + 1); // F1*G1
	add_poly_to(tmp[0], tmp[1], 0, k);
	add_poly_to(tmp[2], tmp[3], 0, k);
	caratsuba_multiplication(tmp[0], tmp[2], tmp[1], k, lvl + 1); // (F0 + F1)*(G0 + G1)
	std::fill(tmp[2].begin(), tmp[2].begin() + (length << 1), 0);
	add_poly_to(tmp[2], tmp[5], length, length);
	add_poly_to(tmp[2], tmp[5], k, length);
	add_poly_to(tmp[2], tmp[4], k, length);
	add_poly_to(tmp[2], tmp[1], k, length);
	add_poly_to(tmp[2], tmp[4], 0, length);
	std::copy(tmp[2].begin(), tmp[2].begin() + std::min((size_t)length << 1, dst.size()), dst.begin());
	return dst;
}

std::vector<unsigned>& galois_field::SCHONHAGE_STRASSEN_FFT(std::vector<unsigned>& a, std::vector<unsigned>& b, std::vector<unsigned>& dst, unsigned n, unsigned level) {
	++_poly_multiplications;
	unsigned c = 1, k = 0;
	while (c < n) {
		c *= 3;
		++k;
	}
	size_t start_mult = _multiplications;
	size_t start_add = _additions;
	unsigned m = std::pow(3, k / 2 + k % 2);
	auto& tmp = _schonhage_strassen_tmp[level];
	for (auto& v : tmp) {
		std::fill(v.begin(), v.end(), 0);
	}

	std::copy(a.begin(), a.begin() + std::min((size_t)2 * n, a.size()), tmp[6].begin());
	std::copy(b.begin(), b.begin() + std::min((size_t)2 * n, b.size()), tmp[7].begin());

	if (n <= 27) {
		for (size_t i = 0; i < 2 * n; ++i) {
			for (size_t j = 0; j < 2 * n; ++j) {
				if (i + j >= dst.size()) {
					break;
				}
				tmp[2][i + j] = add(tmp[2][i + j], multiply(tmp[6][i], tmp[7][j]));
			}
		}
		add_subpoly_with_modular_shift(tmp[2], tmp[2], 0, 2 * n, 4 * n, n, 2 * n);
		std::copy(tmp[2].begin(), tmp[2].begin() + std::min((size_t)2 * n, dst.size()), dst.begin());
		return dst;
	}
	if (n <= 256) {
		size_t len = 1 << (sizeof(unsigned) * 8 - std::countl_zero(n - 1) + 1);
		caratsuba_multiplication(tmp[6], tmp[7], tmp[2], len, 0);
		add_subpoly_with_modular_shift(tmp[2], tmp[2], 0, 2 * n, 4 * n, n, 2 * n);
		std::copy(tmp[2].begin(), tmp[2].begin() + std::min((size_t)2 * n, dst.size()), dst.begin());
		return dst;
	}
	unsigned t = n / m;
	unsigned eta = 1;
	unsigned block_size = 2 * m;
	if (t != m) {
		eta = 3;
	}


	for (size_t i = 1; i <= 2; ++i) {

		for (size_t _ = 0; _ < 4; ++_) {
			std::fill(tmp[_].begin(), tmp[_].end(), 0);
		}
		// preparations
		for (size_t j = 0; j < t; ++j) {
			// f*_i <- f* - eta^(it)
			std::copy(tmp[6].begin() + block_size * j / 2, tmp[6].begin() + block_size * (j + 1) / 2, tmp[0].begin() + block_size * j);
			add_subpoly_with_modular_shift(tmp[0], tmp[6], block_size * j, n + block_size * j / 2, n + block_size * (j + 1) / 2, m, eta * i * t);


			// g*_i <- g* - eta^(it)
			std::copy(tmp[7].begin() + block_size * j / 2, tmp[7].begin() + block_size * (j + 1) / 2, tmp[1].begin() + block_size * j);
			add_subpoly_with_modular_shift(tmp[1], tmp[7], block_size * j, n + block_size * j / 2, n + block_size * (j + 1) / 2, m, eta * i * t);
			// substitute y -> (eta^i)*y
			add_subpoly_with_modular_shift(tmp[2], tmp[0], block_size * j, block_size * j, block_size * (j + 1), m, eta * i * j);
			add_subpoly_with_modular_shift(tmp[3], tmp[1], block_size * j, block_size * j, block_size * (j + 1), m, eta * i * j);

		}

		// call convolution
		SCHONHAGE_CONVOLUTION(tmp[2], tmp[3], tmp[3 + i], t, eta * 3, m, level);
	}

	std::fill(tmp[0].begin(), tmp[0].end(), 0);
	std::fill(tmp[1].begin(), tmp[1].end(), 0);

	for (size_t j = 0; j < t; ++j) {
		// backwards substitution y -> (eta^-i)*y
		add_subpoly_with_modular_shift(tmp[0], tmp[4], block_size * j, block_size * j, block_size * (j + 1), m, 3 * m - (eta * j) % (3 * m));
		add_subpoly_with_modular_shift(tmp[1], tmp[5], block_size * j, block_size * j, block_size * (j + 1), m, 3 * m - (2 * eta * j) % (3 * m));

	}
	using std::swap;
	swap(tmp[0], tmp[4]);
	swap(tmp[1], tmp[5]);
	std::fill(tmp[0].begin(), tmp[0].end(), 0);
	std::fill(tmp[1].begin(), tmp[1].end(), 0);

	// restore the answer by formula
	// 4 5

	for (size_t j = 0; j < t; ++j) {
		// ?
		add_subpoly_with_modular_shift(tmp[0], tmp[4], 2 * n + block_size * j, block_size * j, block_size * (j + 1), m, 0);
		add_subpoly_with_modular_shift(tmp[0], tmp[5], 2 * n + block_size * j, block_size * j, block_size * (j + 1), m, 0);
		add_subpoly_with_modular_shift(tmp[0], tmp[4], block_size * j, block_size * j, block_size * (j + 1), m, 2 * eta * t);
		add_subpoly_with_modular_shift(tmp[0], tmp[5], block_size * j, block_size * j, block_size * (j + 1), m, eta * t); 
	}

	std::fill(tmp[1].begin(), tmp[1].end(), 0);
	// substitute y = x^m

	for (size_t i = 0; i < 2 * t; ++i) {
		for (size_t j = 0; j < block_size; ++j) {
			size_t id = i * m + j;
			tmp[1][id] = add(tmp[1][id], tmp[0][i * block_size +j]);
		}
	}
	// in modulo x^2n + x^n + 1
	add_subpoly_with_modular_shift(tmp[1], tmp[1], 0, 2 * n, 4 * n, n, 2 * n);
	std::copy(tmp[1].begin(), tmp[1].begin() + std::min((size_t)2ull * n, dst.size()), dst.begin());
	return dst;
}

// basis transform perform is needed, we assume src is already in the correct basis
// such dft allows to decrease overall decoding complexity to nlog(n)log(log(n)) for binary galois field (non-fft-friendly field case)
std::vector<unsigned>& galois_field::DFTimpl(std::vector<unsigned>& src, std::vector<unsigned>& dst, unsigned size, unsigned i, unsigned r) {
	unsigned k = 8 * sizeof(unsigned) - std::countl_zero<unsigned>(size) - 1;
	if (i == k) {
		dst[0] = src[r]; // !!
		return dst;
	}
	// call dfts
	/*std::fill(_dft_tmp[i + 1][r + (1 << i)].begin(), _dft_tmp[i + 1][r + (1 << i)].end(), 0);
	std::fill(_dft_tmp[i + 1][r].begin(), _dft_tmp[i + 1][r].end(), 0);*/
	DFTimpl(src, _dft_tmp[i + 1][r], size, i + 1, r);
	DFTimpl(src, _dft_tmp[i + 1][r + (1 << i)], size, i + 1, r + (1 << i));
	for (unsigned j = 0; j < (1 << (k - i - 1)); ++j) {
		unsigned pos = j << (i + 1);
		auto a = _dft_tmp[i + 1][r][pos];
		auto b = _dft_tmp[i + 1][r + (1 << i)][pos];
		dst[pos] = add(a, multiply(_s[i][pos], b));
		dst[pos + (1 << i)] = add(dst[pos], b);
	}
	return dst;
} // binary architecture a with variable size


std::vector<unsigned>& galois_field::IDFTimpl(std::vector<unsigned>& src, std::vector<unsigned>& dst, unsigned size, unsigned i, unsigned r) {
	unsigned k = 8 * sizeof(unsigned) - std::countl_zero<unsigned>(size) - 1;
	if (i == k) {
		dst[r] = src[0];
		return dst;
	}
	std::fill(_dft_tmp[i + 1][r + (1 << i)].begin(), _dft_tmp[i + 1][r + (1 << i)].end(), 0);
	std::fill(_dft_tmp[i + 1][r].begin(), _dft_tmp[i + 1][r].end(), 0);
	/*std::fill(_idft_tmp[i + 1][r + (1 << i)].begin(), _idft_tmp[i + 1][r + (1 << i)].end(), 0);
	std::fill(_idft_tmp[i + 1][r].begin(), _idft_tmp[i + 1][r].end(), 0);*/
	for (unsigned j = 0; j < (1 << (k - i - 1)); ++j) {
		unsigned pos = j << (i + 1);
		_dft_tmp[i + 1][r + (1 << i)][pos] = add(src[pos], src[pos + (1 << i)]);
		_dft_tmp[i + 1][r][pos] = add(src[pos], multiply(_s[i][pos], _dft_tmp[i + 1][r + (1 << i)][pos]));
	}

	IDFTimpl(_dft_tmp[i + 1][r], _idft_tmp[i + 1][r], size, i + 1, r);
	IDFTimpl(_dft_tmp[i + 1][r + (1 << i)], _idft_tmp[i + 1][r + (1 << i)], size, i + 1, r + (1 << i));
	// count an answer
	/*for (unsigned j = 0; j < (1 << (k - i - 1)); ++j) {
		dst[2 * j] = _idft_tmp[i + 1][r][j];
		dst[2 * j + 1] = _idft_tmp[i + 1][r + (1 << i)][j];

	}*/
	for (unsigned j = 0; j < (1 << (k - i - 1)); ++j) {
		unsigned pos = (j << (i + 1)) + r;
		dst[pos] = _idft_tmp[i + 1][r][pos];
		dst[pos + (1 << i)] = _idft_tmp[i + 1][r + (1 << i)][(pos + (1 << i))];	
	}
	return dst;
} // binary architecture a with variable size

void galois_field::split_poly(std::vector<unsigned> const& p, std::vector<unsigned>& dst1, std::vector<unsigned>& dst2, unsigned m) {
	std::fill(dst2.begin(), dst2.end(), 0);
	std::fill(dst1.begin(), dst1.end(), 0);
	for (unsigned i = 0; i <= degree(p); ++i) {
		if (i < m) {
			dst2[i] = p[i];
		}
		else
		{
			dst1[i - m] = p[i];
		}
	}
}

std::vector<unsigned>& galois_field::add_poly(std::vector<unsigned>& a, std::vector<unsigned>& b, std::vector<unsigned>& dst, unsigned m) {
	std::fill(dst.begin(), dst.end(), 0);
	for (unsigned i = 0; i < a.size(); ++i) {
		if (i < m) {
			dst[i] = a[i];
		}
		else {
			dst[i] = add(a[i], b[i - m]);
		}

	}
	return dst;
}

std::vector<unsigned>& galois_field::add_poly_to(std::vector<unsigned>& a, std::vector<unsigned>& b, unsigned m, unsigned length) {
	for (unsigned i = 0; i < std::min(b.size(), (size_t)length); ++i) {
		if (i + m >= a.size()) {
			throw std::runtime_error("polynomial overflow");
		}
		a[i + m] = add(a[i + m], b[i]);

	}
	return a;
}

std::vector<unsigned>& galois_field::sub_poly(std::vector<unsigned>& a, std::vector<unsigned>& b, std::vector<unsigned>& dst) {
	unsigned i;
	for (i = 0; i < a.size() && i < b.size(); ++i) {
		dst[i] = add(a[i], b[i]);
	}
	while (i < a.size()) {
		dst[i] = a[i];
		++i;
	}
	while (i < b.size()) {
		dst[i] = b[i];
		++i;
	}
	return dst;
}

std::vector<unsigned>& galois_field::multipy_poly_by_const(std::vector<unsigned>& p, unsigned a) {
	for (unsigned i = 0; i < p.size(); ++i) {
		p[i] = multiply(p[i], a);
	}
	return p;
}

std::vector<unsigned>& galois_field::rev_poly(std::vector<unsigned>& src, std::vector<unsigned>& dst, unsigned n) {
	for (unsigned i = 0; i < n + 1; ++i) {
		dst[n - i] = src[i];
	}
	std::fill(dst.begin() + n + 1, dst.end(), 0);
	return dst;
}

void galois_field::print_poly(std::vector<unsigned> const& p) {
	for (auto v : p) {
		std::cout << v << " ";
	}
	std::cout << "\n";

}

std::vector<unsigned>& galois_field::inv_poly(std::vector<unsigned>& src, std::vector<unsigned>& dst, unsigned mod) {
	// make zeros in inverse tmp
	unsigned g0 = inverse(src[0]);
	unsigned r = sizeof(unsigned) * 8 - std::countl_zero(mod - 1);
	dst[0] = g0;
	size_t start_add = _additions;
	size_t start_mult = _multiplications;
	for (unsigned i = 1; i <= r; ++i) {
		fast_poly_multiplication(dst, dst, _inverse_temporary1);
		remainder_of_power(_inverse_temporary1, 1 << (i));
		std::copy(src.begin(), src.begin() + (1 << i), _inverse_temporary2.begin());

		remainder_of_power(_inverse_temporary2, 1 << (i));
		fast_poly_multiplication(_inverse_temporary2, _inverse_temporary1, dst, (unsigned)(1 << (i)));
		remainder_of_power(dst, std::min((unsigned)(1 << i), mod));
	}
	return dst;
}

std::vector<unsigned>& galois_field::remainder_of_power(std::vector<unsigned>& a, unsigned n) {
	for (unsigned i = n; i < a.size(); ++i) {
		a[i] = 0;
	}
	return a;
}

std::vector<unsigned>& galois_field::formal_derivate(std::vector<unsigned>& src, std::vector<unsigned>& dst) {
	std::fill(dst.begin(), dst.end(), 0);
	for (size_t i = 1; i < src.size(); i += 2) {
		dst[i - 1] = src[i];
	}
	return dst;
}

std::vector<unsigned>& galois_field::translate_bit_vector(linalg::bit_vector& src, std::vector<unsigned>& dst) {
	unsigned cur = 0;
	for (size_t i = 0; i < src.size(); ++i) {
		cur += src[i] << (i % _m);
		if (i % _m == _m - 1) {
			dst[i / _m] = cur;
			cur = 0;
		}
	}
	if (src.size() % _m != 0) {
		dst[src.size() / _m + 1] = cur;
	}
	return dst;

}
linalg::bit_vector& galois_field::translate_to_bit_vector(std::vector<unsigned>& src, linalg::bit_vector& dst) {
	unsigned cur = 0;
	for (size_t i = 0; i < dst.size(); ++i) {
		if (i % _m == 0) {
			cur = src[i / _m];
		}
		dst.set(i, (cur >> (i % _m)) & 1);
	}
	return dst;
}

unsigned galois_field::substitute_poly(std::vector<unsigned>& poly, unsigned val) {
	unsigned result = 0;
	for (ptrdiff_t i = poly.size() - 1; i >= 0; --i) {
		result = multiply(result, val);
		result = add(result, poly[i]);

	}
	return result;

}

// assume t always equals 2, n is always a power of 2
std::vector<unsigned>& galois_field::taylor_expansion(std::vector<unsigned>& src, std::vector<unsigned>& dst, unsigned n, unsigned t) {
	if (n <= t) {
		std::copy_n(src.begin(), n, dst.begin());
		return dst;
	}
	std:fill(dst.begin(), dst.begin() + n, 0);
	unsigned k = sizeof(unsigned) * 8 - std::countl_zero((n / t) - 1) - 1;
	auto& tmp = _taylor_expansion_tmp[k + 2];
	for (auto& v : tmp) {
		std::fill(v.begin(), v.begin() + n, 0);
	}
	unsigned p = t << k;
	unsigned r = (t - 1) << k;
	// compute h
	if (n - p < r) {
		add_subpoly(tmp[0], src, 0, p, n); // tmp[0] := h
	}
	else {
		add_subpoly(tmp[0], src, 0, p, p + r); // tmp[0] := h
		add_subpoly(tmp[0], src, 0, p + r, n); 
	}
	// compute g0
	add_subpoly(tmp[1], src, 0, 0, p);
	add_subpoly(tmp[1], tmp[0], 1 << k, 0, r);

	// first part of answer
	taylor_expansion(tmp[1], tmp[2], p, t);

	// compute g1
	std::fill(tmp[1].begin(), tmp[1].end(), 0);
	add_subpoly(tmp[1], tmp[0], 0, 0, r);
	add_subpoly(tmp[1], src, r, p + r, n); // <- error is probable here
	std::fill(tmp[0].begin(), tmp[0].end(), 0);
	taylor_expansion(tmp[1], tmp[0], n - p, t);
	// (g0 + g1x) + (g2 + g3x)(x + x^2) + (g4 + g5x)(x + x^2)^2..
	// do result
	std::copy(tmp[2].begin(), tmp[2].begin() + p, dst.begin());
	std::copy(tmp[0].begin(), tmp[0].begin() + n - p, dst.begin() + p);
	return dst;
}

// assume t always equals 2, n is always a power of 2
// it is reversed order taylor_expansion
std::vector<unsigned>& galois_field::itaylor_expansion(std::vector<unsigned>& src, std::vector<unsigned>& dst, unsigned n, unsigned t) {
	if (n <= t) {
		std::copy_n(src.begin(), n, dst.begin());
		return dst;
	}
	std:fill(dst.begin(), dst.begin() + n, 0);
	unsigned k = sizeof(unsigned) * 8 - std::countl_zero((n / t) - 1) - 1;
	auto& tmp = _taylor_expansion_tmp[k + 2];
	for (auto& v : tmp) {
		std::fill(v.begin(), v.begin() + n, 0);
	}
	unsigned p = t << k;
	unsigned r = (t - 1) << k;
	std::copy(src.begin(), src.begin() + p, tmp[0].begin());
	itaylor_expansion(tmp[0], tmp[2], p, t); // g0
	std::copy(src.begin() + p, src.begin() + n, tmp[1].begin());
	itaylor_expansion(tmp[1], tmp[0], n - p, t); // g1
	add_subpoly(tmp[2], tmp[0], 1 << k, 0, r);
	add_subpoly(tmp[0], tmp[0], 0, r, n - p);
	std::copy(tmp[2].begin(), tmp[2].begin() + p, dst.begin());
	std::copy(tmp[0].begin(), tmp[0].begin() + n - p, dst.begin() + p);
	return dst;
}


std::vector<unsigned>& galois_field::gao_mateer_fft(std::vector<unsigned>& src, std::vector<unsigned>& dst, unsigned m) {
	if (m == 1) {
		dst[0] = src[0];
		dst[1] = add(src[0], multiply(src[1], precomputed_basises_delta[m][0]));

		return dst;
	}
	unsigned n = 1 << m;
	auto& tmp = _gao_mateer_fft_tmp[m];
	// substitute x=beta_m*x
	unsigned beta = 1;
	for (size_t i = 0; i < n; ++i) {
		tmp[0][i] = multiply(src[i], beta);
		beta = multiply(beta, precomputed_basises_delta[m][m - 1]);
	}
	// taylor expansion + split
	taylor_expansion(tmp[0], tmp[1], n, 2);
	for (size_t i = 0; i < n / 2; ++i) {
		tmp[2][i] = tmp[1][2 * i];
		tmp[3][i] = tmp[1][2 * i + 1];
	}
	// pre answers
	gao_mateer_fft(tmp[2], tmp[0], m - 1);
	gao_mateer_fft(tmp[3], tmp[1], m - 1);
	// combine pre answers
	// G[i] = alpha_0*gamma_0 + .. + alpha_m-1*gamma_m-1 = (alpha_0beta_0 + .. + alpha_m-1beta_m-1)beta_m^-1
	unsigned k = 1 << (m - 1);
	for (size_t i = 0; i < k; ++i) {
		dst[i] = add(tmp[0][i], multiply(precomputed_space_gamma[m - 1][i], tmp[1][i]));
		dst[k + i] = add(dst[i], tmp[1][i]);
	}
	return dst;
}

std::vector<unsigned>& galois_field::gao_mateer_ifft(std::vector<unsigned>& src, std::vector<unsigned>& dst, unsigned m) {
	if (m == 1) {
		dst[0] = src[0];
		dst[1] = multiply(add(src[0], src[1]), inverse(precomputed_basises_delta[m][0]));
		return dst;
	}
	unsigned n = 1 << m;
	unsigned k = 1 << (m - 1);
	auto& tmp = _gao_mateer_fft_tmp[m];
	for (size_t i = 0; i < k; ++i) {
		tmp[1][i] = add(src[i], src[k + i]);
		tmp[0][i] = (i == 0 ? src[i] : add(src[i], multiply(inverse(precomputed_space_gamma[m - 1][i]), tmp[1][i])));
	}

	gao_mateer_ifft(tmp[0], tmp[2], m - 1);
	gao_mateer_ifft(tmp[1], tmp[3], m - 1);

	for (size_t i = 0; i < n / 2; ++i) {
		tmp[1][2 * i] = tmp[2][i];
		tmp[1][2 * i + 1] = tmp[3][i];
	}
	
	itaylor_expansion(tmp[1], tmp[0], n, 2);
	print_poly(tmp[1]);
	print_poly(tmp[0]);
	unsigned beta = 1;
	unsigned step = inverse(precomputed_basises_delta[m][m - 1]);
	for (size_t i = 0; i < n; ++i) {
		dst[i] = multiply(tmp[0][i], beta);
		beta = multiply(beta, step);
	}
	
	return dst;
}
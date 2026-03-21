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
	, _division_tmp(_n + 5)
	, _const2(_n + 5, 0)
	, _a_tmp(_q)
	, _b_tmp(_q)
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

	// evaluate s
	// recurrent way of evaluation may be easier
	_s.resize(_m + 1);
	_s.back().resize(_q);
	_s[0].resize(_q);
	for (size_t i = 1; i < _q; ++i) {
		_s[0][i] = i;
	}
	for (size_t level = 1; level < _m; ++level) {
		_s[level].resize(_q);
		for (size_t i = 0; i < _q; ++i) {
			_s[level][i] = 1;
			for (size_t k = 0; k < (1 << level); ++k) {
				_s[level][i] = multiply(_s[level][i], add(k, i));
			}
		}
	}
	std::cout << "s:" << _s.size() << "\n";
	for (size_t level = 0; level < _m; ++level) {
		std::cout << level << ": ";
		unsigned inv = inverse(_s[level][1 << level]);
		for (size_t i = 0; i < _q; ++i) {
			_s[level][i] = multiply(_s[level][i], inv);
			std::cout << _s[level][i] << " ";
		}
		std::cout << "\n";
	}

	std::cout << 3 << ": ";
	for (size_t i = 0; i < _q; ++i) {
		std::cout << _s[3][i] << " ";
	}
	std::cout << "\n";
	size_t tmp_sizes = _q;
	_emgcd_tmp_polynomials.resize(40);
	for (auto& tmp : _emgcd_tmp_polynomials) {
		for (auto& v : tmp) {
			v.resize(tmp_sizes);
		}
	}

	_emgcd_tmp_result.resize(40);
	for (auto& tmp : _emgcd_tmp_result) {
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
	for (auto& tmp : _solve_toeplitz_tmp) {
		//for (auto& v : tmp) {
		tmp.resize(tmp_sizes);
			//v.second.resize(_n);
		//}
	}

	std::cout << "initted\n";
}

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
	for (ptrdiff_t deg = a.size() - 1; deg >= 0; --deg) {
		if (a[deg] != 0) {
			return deg;
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
	return a ^ b;
}

unsigned galois_field::multiply(unsigned a, unsigned b) const { 
	if (a == 0 || b == 0) {
		return 0;
	}
	unsigned log_sum = (_log_table[a] + _log_table[b]) % _n;
	return _exp_table[log_sum]; 
}
unsigned galois_field::inverse(unsigned a) const {
	assert(a != 0);
	if (a == 0) throw std::runtime_error("division by zero");
    return _exp_table[_n - _log_table[a]];
}
unsigned galois_field::divide(unsigned a, unsigned b) const { 
	return multiply(a, inverse(b)); 
}

std::vector<unsigned>& galois_field::fast_poly_multiplication(std::vector<unsigned>& a, std::vector<unsigned>& b, std::vector<unsigned>& dst) {
	/*DFT(a, _a_tmp);
	DFT(b, _b_tmp);
	
	for (unsigned i = 0; i < a.size(); ++i) {
		_a_tmp[i] = multiply(_a_tmp[i], _b_tmp[i]);
	}
	
	IDFT(_a_tmp, dst);*/
	for (auto& v : dst) {
		v = 0;
	}
	for (size_t i = 0; i < a.size(); ++i) {
		for (size_t j = 0; j < b.size(); ++j) {
			if (i + j >= dst.size()) {
				break;
			}
			dst[i + j] = add(dst[i + j], multiply(a[i], b[j]));
		}
	}
	//std::reverse(dst.begin() + 1, dst.end());
	//multipy_poly_by_const(dst, _inverse_element);

	return dst;
}

std::vector<unsigned>& galois_field::fast_poly_division(std::vector<unsigned>& a, std::vector<unsigned>& b, std::vector<unsigned>& quotient, std::vector<unsigned>& remainder) {
	auto deg_a = degree(a), deg_b = degree(b);
	std::fill(remainder.begin(), remainder.end(), 0);
	std::fill(quotient.begin(), quotient.end(), 0);
	if (deg_a < deg_b) {
		std::copy(a.begin(), a.end(), remainder.begin());

		return a;
	}
	auto m = deg_a - deg_b;
	//std::cout << deg_a << " " << deg_b << "\n";
	//std::cout << "ARGS:\n";
	//print_poly(a);
	//print_poly(b);
	//print_poly(_division_tmp);
	//std::cout << "\n";

	rev_poly(a, quotient, deg_a);
	rev_poly(b, _division_tmp, deg_b);
	//print_poly(_division_tmp);
	inv_poly(_division_tmp, remainder, m + 1);
	//print_poly(remainder);
	//std::cout << "here\n";
	fast_poly_multiplication(quotient, remainder, _division_tmp);
	//print_poly(_division_tmp);
	remainder_of_power(_division_tmp, m + 1);
	//print_poly(_division_tmp);
	rev_poly(_division_tmp, quotient, m);
	fast_poly_multiplication(quotient, b, _division_tmp);
	//print_poly(quotient);
	//print_poly(_division_tmp);
	add_poly(a, _division_tmp, remainder, 0);
	return quotient;
}

std::array<std::pair<std::vector<unsigned>, std::vector<unsigned>>, 3>& 
		galois_field::EMGCD(std::vector<unsigned> const& u0, std::vector<unsigned> const& u1, std::array<std::pair<std::vector<unsigned>, std::vector<unsigned>>, 3>& dst, unsigned tmp_num) {
	unsigned n = degree(u0);
	//std::cout << "enetring emgcd with args:\n";
	//print_poly(u0);
	//print_poly(u1);
	auto deg_u1 = degree(u1);
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
		return dst;
	}
	auto& locals = _emgcd_tmp_polynomials[tmp_num];
	for (auto& v : locals) {
		std::fill(v.begin(), v.end(), 0);
	}
	std::array<std::pair<std::vector<unsigned>, std::vector<unsigned>>, 3>& result = _emgcd_tmp_result[tmp_num];
	//std::vector<unsigned> b0, b1, c0, c1;
	// tmp layout:
	// b0   b1   c0   c1   d    e    q    f    
	// 0    1    2    3    4    5    6    7  
	//std::cout << tmp_num << "\n";
	//print_poly(u0);
	//print_poly(u1);
	split_poly(u0,locals[0], locals[2], m);
	split_poly(u1, locals[1], locals[3], m);
	//print_poly(locals[0]);
	//std::cout << "here "  << tmp_num << " " << n << " " << m << " " << deg_u1 << "\n";
	// u1w1v1 is stored in _emgcd_tmp_result
	EMGCD(locals[0], locals[1], _emgcd_tmp_result[tmp_num], tmp_num + 1);

	using std::swap;
	//std::vector<unsigned> d,e;
	//
	
	//std::cout << "here1\n";
	// b0, b1 wont be used further
	fast_poly_multiplication(result[1].first, locals[2], locals[4]);
	fast_poly_multiplication(result[2].first, locals[3], locals[1]);
	
	add_poly(locals[4], locals[1], locals[0], 0);
	//swap(locals[0], locals[4]);
	add_poly(locals[0], result[0].first, locals[4], m); // d
	//print_poly(locals[4]);
	// in some locals may be litter
	fast_poly_multiplication(result[1].second, locals[2], locals[5]);
	fast_poly_multiplication(result[2].second, locals[3], locals[1]);
	add_poly(locals[5], locals[1], locals[0], 0);
	//swap(locals[0], locals[5]);
	add_poly(locals[0], result[0].second, locals[5], m); // e
	//print_poly(locals[5]);

	unsigned deg_e = degree(locals[5]);
	//std::cout << deg_e << " " << m << "\n";
	if (deg_e < m) {
		dst[0].first = std::move(locals[0]);
		dst[0].second = std::move(locals[1]);
		dst[1].first = std::move(result[1].first);
		dst[1].second = std::move(result[1].second);
		dst[2].first = std::move(result[2].first);
		dst[2].second = std::move(result[2].second);
		return dst;
	}
	//std::cout << "here\n";

	// c0 and c1 also can be reused -> g0,g1,h0,h1 will be stored there (in b0,b1,c0,c1)

	// tmp layout:
	// g0   g1   h0   h1   d    e    q    f    
	// 0    1    2    3    4    5    6    7  
	fast_poly_division(locals[4], locals[5], locals[6], locals[7]);
	unsigned k = 2 * m - deg_e;
	split_poly(locals[5], locals[0], locals[2], k);
	split_poly(locals[7], locals[1], locals[3], k);
	EMGCD(locals[0], locals[1], _emgcd_tmp_result[tmp_num + 1], tmp_num + 1); // need to be careful with such storage

	std::array<std::pair<std::vector<unsigned>, std::vector<unsigned>>, 3>& result2 = _emgcd_tmp_result[tmp_num + 1];
	
	// g0, g1 wont be used further
	
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
	fast_poly_multiplication(result[2].first, locals[6], locals[0]);
	add_poly(result[1].first, locals[0], locals[1], 0);

	fast_poly_multiplication(result[2].second, locals[6], locals[0]);
	add_poly(result[1].second, locals[0], locals[2], 0);

	fast_poly_multiplication(result[2].first, result2[1].first, locals[0]);
	fast_poly_multiplication(locals[1], result2[1].second, locals[3]);
	add_poly(locals[0], locals[3], dst[1].first, 0);

	fast_poly_multiplication(result[2].second, result2[1].first, locals[0]);
	fast_poly_multiplication(locals[2], result2[1].second, locals[3]);
	add_poly(locals[0], locals[3], dst[2].first, 0);

	fast_poly_multiplication(result[2].first, result2[2].first, locals[0]);
	fast_poly_multiplication(locals[1], result2[2].second, locals[3]);
	add_poly(locals[0], locals[3], dst[1].second, 0);

	fast_poly_multiplication(result[2].second, result2[2].first, locals[0]);
	fast_poly_multiplication(locals[2], result2[2].second, locals[3]);
	add_poly(locals[0], locals[3], dst[2].second, 0);
	// end multiply 


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
	//std::cout << "here?\n";
	_ad_tmp_polinomyals[0][2 * n + 1] = 0;
	_ad_tmp_polinomyals[0][2 * n + 3] = 1;
	shift_poly(_ad_tmp_polinomyals[1]);
	_ad_tmp_polinomyals[1][2 * n + 2] = 1;
	_ad_tmp_polinomyals[1][0] = 1;
	EMGCD(_ad_tmp_polinomyals[0], _ad_tmp_polinomyals[1], _ad_tmp_emgcd_results[0], 0);
	//std::cout << degree(_ad_tmp_emgcd_results[0][0].second) << " " << n + 1 << "\n";
	if (degree(_ad_tmp_emgcd_results[0][0].second) == n + 1) {
		rev_poly(_ad_tmp_polinomyals[1], _ad_tmp_polinomyals[2], 2 * n + 2);
		EMGCD(_ad_tmp_polinomyals[0], _ad_tmp_polinomyals[2], _ad_tmp_emgcd_results[1], 0);
		dst1 = multipy_poly_by_const(_ad_tmp_emgcd_results[0][2].second, inverse(_ad_tmp_emgcd_results[0][0].second[degree(_ad_tmp_emgcd_results[0][0].second)]));
		dst2 = multipy_poly_by_const(_ad_tmp_emgcd_results[1][2].second, inverse(_ad_tmp_emgcd_results[1][0].second[degree(_ad_tmp_emgcd_results[1][0].second)]));
		return;
	}
	//std::cout << "hhhhh!\n";
	//assert(false);
	_ad_tmp_polinomyals[1][2 * n + 2] = _inverse_element;
	//print_poly(_ad_tmp_polinomyals[0]);
	//print_poly(_ad_tmp_polinomyals[1]);
	EMGCD(_ad_tmp_polinomyals[0], _ad_tmp_polinomyals[1], _ad_tmp_emgcd_results[0], 0);
	rev_poly(_ad_tmp_polinomyals[1], _ad_tmp_polinomyals[2], 2 * n + 2);
	//print_poly(_ad_tmp_polinomyals[0]);
	//print_poly(_ad_tmp_polinomyals[2]);
	EMGCD(_ad_tmp_polinomyals[0], _ad_tmp_polinomyals[2], _ad_tmp_emgcd_results[1], 0);
	//std::cout << "here&&&\n";
	dst1 = multipy_poly_by_const(_ad_tmp_emgcd_results[0][2].second, inverse(_ad_tmp_emgcd_results[0][0].second[degree(_ad_tmp_emgcd_results[0][0].second)]));
	dst2 = multipy_poly_by_const(_ad_tmp_emgcd_results[1][2].second, inverse(_ad_tmp_emgcd_results[1][0].second[degree(_ad_tmp_emgcd_results[1][0].second)]));
}

std::vector<unsigned>& galois_field::SOLVE_TOEPITZ(std::vector<unsigned>& a, std::vector<unsigned>& b, unsigned n, std::vector<unsigned>& dst) {
	AD(a, n, _ad_x, _ad_y);
	for (auto& v : _solve_toeplitz_tmp) {
		std::fill(v.begin(), v.end(), 0);
	}
	//std::cout << "here__\n";
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
	call_fft(src, dst);
	return dst;
}

std::vector<unsigned>& galois_field::IDFT(std::vector<unsigned>& src, std::vector<unsigned>& dst) {
	DFT(src, dst);
	/*for (unsigned i = 0; i < dst.size(); ++i) {
		dst[i] = multiply(dst[i], _inverse_element);
	}*/
	std::reverse(dst.begin() + 1, dst.end() - 1);
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

std::vector<unsigned>& galois_field::SCHONHAGE_STRASSEN_FFT(std::vector<unsigned>& src, std::vector<unsigned>& dst, unsigned n) {
	if (n <= 2) {

	}
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
	for (unsigned i = 0; i < p.size(); ++i) {
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
	//std::cout << "add\n";
	//print_poly(dst);
	return dst;
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

std::vector<unsigned> galois_field::inv_poly(std::vector<unsigned>& src, std::vector<unsigned>& dst, unsigned mod) {
	// make zeros in inverse tmp
	unsigned g0 = inverse(src[0]);
	unsigned r = sizeof(unsigned) * 8 - std::countl_zero(mod - 1);
	dst[0] = g0;
	for (unsigned i = 1; i <= r; ++i) {
		fast_poly_multiplication(dst, dst, _inverse_temporary1);
		remainder_of_power(fast_poly_multiplication(src, _inverse_temporary1, dst), std::min((unsigned)(1 << r), mod));
	}

	return dst;
}

std::vector<unsigned>& galois_field::remainder_of_power(std::vector<unsigned>& a, unsigned n) {
	for (unsigned i = n; i < a.size(); ++i) {
		a[i] = 0;
	}
	return a;
}


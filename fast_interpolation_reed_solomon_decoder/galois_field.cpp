#include "galois_field.h"
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <vector>
#include <bit>
#include <cassert>

#include <iostream>

galois_field::galois_field(size_t p, size_t m, std::vector<size_t> const& gen) 
	: _p(p)
	, _m(m)
	, _q(std::pow(p, m))
	, _n(_q - 1)
	, _generating_poly(gen)
	, _log_table(_q) 
	, _exp_table(_q)
	, _minus_one(_p - 1)
	, _a_tmp(_n)
	, _b_tmp(_n)
{
	init();
}

void galois_field::init() {
	std::vector<size_t> x(_m + 1);
	x[0] = 1;
	for (size_t i = 0; i < _q; ++i) {
		_exp_table[i] = poly_to_num(x);
		_log_table[poly_to_num(x)] = i;
		shift_poly(x);
		if (x.back()) {
			x = add_init(x, inverse_add_init(multiply_by_const_init(_generating_poly, x.back())));
		}
	}
	_exp_table[_n] = _exp_table[0];
}

std::vector<size_t> galois_field::inverse_add_init(std::vector<size_t> const& v) const {
	std::vector<size_t> res;
	for (size_t i = 0; i < v.size(); ++i) {
		res.emplace_back((_p - v[i]) % _p);
	}
	return  res;
}

std::vector<size_t> galois_field::add_init(std::vector<size_t> const& a, std::vector<size_t> const& b) const {
	std::vector<size_t> res;
	size_t i = 0;
	while (i < a.size() && i < b.size()) {
		res.emplace_back((a[i] + b[i]) % _p);
		++i;
	}
	while (i < a.size()) {
		res.emplace_back(a[i]);
		++i;
	}
	while (i < b.size()) {
		res.emplace_back(b[i]);
		++i;
	}
	return res;
}

std::vector<size_t> galois_field::multiply_by_const_init(std::vector<size_t> const& a, size_t b) const {
	std::vector<size_t> res;
	for (size_t i : a) {
		res.push_back(i * b % _p);
	}
	return res;
}

size_t galois_field::degree(std::vector<size_t> const& a) {
	for (ptrdiff_t deg = a.size() - 1; deg >= 0; --deg) {
		if (a[deg] != 0) {
			return deg;
		}
	}
	return 0;
}

size_t galois_field::poly_to_num(std::vector<size_t> const& a) const {
	size_t res = 0;
	for (ptrdiff_t i = a.size() - 1; i >= 0; --i) {
		res *= _p;
		res += a[i];
	}
	return res;
}

std::vector<size_t>& galois_field::shift_poly(std::vector<size_t>& a) const {
	for (ptrdiff_t i = a.size() - 1; i >= 1; --i) {
		a[i] = a[i - 1];
	}
	a[0] = 0;
	return a;
}

size_t galois_field::add(size_t a, size_t b) const { 
	if (a == 0) return b;
	if (b == 0) return a;
	size_t m, n;
	m = _log_table[a];
	n = _log_table[b];
	std::cout << m << " " << n << "\n";
	if (m <= n) {
		n = _log_table[_exp_table[n - m] + 1];
		size_t log_sum = (m + n) % _n;
		return _exp_table[log_sum];
	}
	else {
		std::cout << _exp_table[m - n] << " " << _log_table[_exp_table[m - n] + 1] << "\n";
		m = _log_table[_exp_table[m - n] + 1];
		size_t log_sum = (m + n) % _n;
		return _exp_table[log_sum];
	}
}
size_t galois_field::sub(size_t a, size_t b) const { 
	return add(a, multiply(b, _minus_one)); 
}

size_t galois_field::multiply(size_t a, size_t b) const { 
	if (a == 0 || b == 0) {
		return 0;
	}
	size_t log_sum = (_log_table[a] + _log_table[b]) % _n;
	return _exp_table[log_sum]; 
}
size_t galois_field::inverse(size_t a) const {
	if (a == 0) throw std::runtime_error("division by zero");
    return _exp_table[_n - _log_table[a]];
}
size_t galois_field::divide(size_t a, size_t b) const { 
	return multiply(a, inverse(b)); 
}

std::vector<size_t>& galois_field::fast_poly_multiplication(std::vector<size_t>& a, std::vector<size_t>& b, std::vector<size_t>& dst) {
	DFT(a, _a_tmp, _exp_table[1]);
	DFT(b, _b_tmp, _exp_table[1]);
	for (size_t i = 0; i < a.size(); ++i) {
		_a_tmp[i] = multiply(_a_tmp[i], _b_tmp[i]);
	}
	IDFT(_a_tmp, dst, inverse(_exp_table[1]));
	return dst;
}

std::vector<size_t>& galois_field::fast_poly_division(std::vector<size_t>& a, std::vector<size_t>& b, std::vector<size_t>& quotient, std::vector<size_t>& remainder) {
	if (a.size() < b.size()) {
		return a;
	}
	//return rev_poly();
}

std::array<std::pair<std::vector<size_t>, std::vector<size_t>>, 3>& 
		galois_field::EMGCD(std::vector<size_t> const& u0, std::vector<size_t> const& u1, std::array<std::pair<std::vector<size_t>, std::vector<size_t>>, 3>& dst, size_t tmp_num) {
	size_t n = degree(u0);

	if (degree(u1) < n / 2) {
		dst[0].first = u0;
		dst[0].second = u1;
		dst[1].first[0] = 1;
		dst[2].second[0] = 1;
		return dst;
	}
	auto& locals = _emgcd_tmp_polynomials[tmp_num];
	std::array<std::pair<std::vector<size_t>, std::vector<size_t>>, 3>& result = _emgcd_tmp_result[tmp_num];
	//std::vector<size_t> b0, b1, c0, c1;
	// tmp layout:
	// b0   b1   c0   c1   d    e    q    f    
	// 0    1    2    3    4    5    6    7  
	size_t m = n / 2 + n % 2;

	split_poly(u0,locals[0], locals[2], m);
	split_poly(u1, locals[1], locals[3], m);

	// u1w1v1 is stored in _emgcd_tmp_result
	EMGCD(locals[0], locals[1], _emgcd_tmp_result[tmp_num], tmp_num + 1);

	using std::swap;
	//std::vector<size_t> d,e;
	//
	
	// b0, b1 wont be used further
	fast_poly_multiplication(result[1].first, locals[2], locals[0]);
	fast_poly_multiplication(result[2].first, locals[3], locals[1]);
	add_poly(locals[0], locals[1], locals[4], 0);
	swap(locals[0], locals[4]);

	add_poly(locals[0], result[0].first, locals[4], m); // d

	// in some locals may be litter
	fast_poly_multiplication(result[1].second, locals[2], locals[0]);
	fast_poly_multiplication(result[2].second, locals[3], locals[1]);
	add_poly(locals[0], locals[1], locals[4], 0);
	swap(locals[0], locals[5]);

	add_poly(locals[0], result[0].second, locals[5], m); // e

	size_t deg_e = degree(locals[5]);
	if (deg_e < n / 2) {
		dst[0].first = std::move(locals[4]);
		dst[0].second = std::move(locals[5]);
		dst[1].first = std::move(result[1].first);
		dst[1].second = std::move(result[1].second);
		dst[2].first = std::move(result[2].first);
		dst[2].second = std::move(result[2].second);
		return dst;
	}

	// c0 and c1 also can be reused -> g0,g1,h0,h1 will be stored there (in b0,b1,c0,c1)

	// tmp layout:
	// g0   g1   h0   h1   d    e    q    f    
	// 0    1    2    3    4    5    6    7  
	fast_poly_division(locals[4], locals[5], locals[6], locals[7]);
	size_t k = 2 * m - deg_e;
	split_poly(locals[5], locals[0], locals[2], k);
	split_poly(locals[7], locals[1], locals[3], k);
	EMGCD(locals[0], locals[1], _emgcd_tmp_result[tmp_num + 1], tmp_num + 1); // need to be careful with such storage

	std::array<std::pair<std::vector<size_t>, std::vector<size_t>>, 3>& result2 = _emgcd_tmp_result[tmp_num + 1];
	
	// g0, g1 wont be used further
	
	// fst part of the answer
	fast_poly_multiplication(result[1].first, locals[2], locals[0]);
	fast_poly_multiplication(result[2].first, locals[3], locals[1]);
	add_poly(locals[0], locals[1], locals[4], 0);
	swap(locals[0], locals[4]);

	add_poly(locals[0], result[0].first, locals[4], k); // d

	// in some locals may be litter
	fast_poly_multiplication(result[1].second, locals[2], locals[0]);
	fast_poly_multiplication(result[2].second, locals[3], locals[1]);
	add_poly(locals[0], locals[1], locals[4], 0);
	swap(locals[0], locals[5]);

	add_poly(locals[0], result[0].second, locals[5], k); // e
	// locals [0..3] can be reused

	// multiply matrices (can be optimized)
	fast_poly_multiplication(result[2].first, locals[6], locals[0]);
	multipy_poly_by_const(locals[0], _p - 1);
	add_poly(result[1].first, locals[0], locals[1], 0);

	fast_poly_multiplication(result[2].second, locals[6], locals[0]);
	multipy_poly_by_const(locals[0], _p - 1);
	add_poly(result[1].second, locals[0], locals[2], 0);

	fast_poly_multiplication(result[2].first, result2[1].first, locals[0]);
	fast_poly_multiplication(locals[1], result2[1].second, locals[3]);
	add_poly(locals[0], locals[3], dst[1].first, 0);

	fast_poly_multiplication(result[2].second, result2[1].first, locals[0]);
	fast_poly_multiplication(locals[2], result2[1].second, locals[3]);
	add_poly(locals[0], locals[3], dst[1].second, 0);

	fast_poly_multiplication(result[2].first, result2[2].first, locals[0]);
	fast_poly_multiplication(locals[1], result2[2].second, locals[3]);
	add_poly(locals[0], locals[3], dst[2].first, 0);

	fast_poly_multiplication(result[2].second, result2[2].first, locals[0]);
	fast_poly_multiplication(locals[2], result2[2].second, locals[3]);
	add_poly(locals[0], locals[3], dst[2].second, 0);
	// end multiply 


	dst[0].first = locals[4];
	dst[0].second = locals[5];

	return dst;
}

void galois_field::AD(std::vector<size_t>& a, size_t n, std::vector<size_t>& dst1, std::vector<size_t>& dst2) {
	_ad_tmp_polinomyals[0][2 * n + 1] = 1;
	_ad_tmp_polinomyals[1] = a;
	EMGCD(_ad_tmp_polinomyals[0], _ad_tmp_polinomyals[1], _ad_tmp_emgcd_results[0], 0);
	if (degree(_ad_tmp_emgcd_results[0][0].second) < n) {
		throw std::runtime_error("singular");
	}
	if (_ad_tmp_emgcd_results[0][2].second[0] != 0) {
		rev_poly(_ad_tmp_polinomyals[1], _ad_tmp_polinomyals[2], 2 * n);
		EMGCD(_ad_tmp_polinomyals[0], _ad_tmp_polinomyals[2], _ad_tmp_emgcd_results[2], 0);
		dst1 = multipy_poly_by_const(_ad_tmp_emgcd_results[0][3].second, inverse(_ad_tmp_emgcd_results[0][0].second[degree(_ad_tmp_emgcd_results[0][0].second)]));
		dst2 = multipy_poly_by_const(_ad_tmp_emgcd_results[1][3].second, inverse(_ad_tmp_emgcd_results[1][0].second[degree(_ad_tmp_emgcd_results[1][0].second)]));
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
		EMGCD(_ad_tmp_polinomyals[0], _ad_tmp_polinomyals[2], _ad_tmp_emgcd_results[2], 0);
		dst1 = multipy_poly_by_const(_ad_tmp_emgcd_results[0][3].second, inverse(_ad_tmp_emgcd_results[0][0].second[degree(_ad_tmp_emgcd_results[0][0].second)]));
		dst2 = multipy_poly_by_const(_ad_tmp_emgcd_results[1][3].second, inverse(_ad_tmp_emgcd_results[1][0].second[degree(_ad_tmp_emgcd_results[1][0].second)]));
		return;
	}

	_ad_tmp_polinomyals[1][2 * n + 2] = _p - 1;
	EMGCD(_ad_tmp_polinomyals[0], _ad_tmp_polinomyals[1], _ad_tmp_emgcd_results[0], 0);
	rev_poly(_ad_tmp_polinomyals[1], _ad_tmp_polinomyals[2], 2 * n + 2);
	EMGCD(_ad_tmp_polinomyals[0], _ad_tmp_polinomyals[2], _ad_tmp_emgcd_results[2], 0);
	dst1 = multipy_poly_by_const(_ad_tmp_emgcd_results[0][3].second, inverse(_ad_tmp_emgcd_results[0][0].second[degree(_ad_tmp_emgcd_results[0][0].second)]));
	dst2 = multipy_poly_by_const(_ad_tmp_emgcd_results[1][3].second, inverse(_ad_tmp_emgcd_results[1][0].second[degree(_ad_tmp_emgcd_results[1][0].second)]));
}

std::vector<size_t>& galois_field::SOLVE_TOEPITZ(std::vector<size_t>& a, std::vector<size_t>& b, size_t n, std::vector<size_t>& dst) {
	AD(a, n, _ad_x, _ad_y);
	assert(_ad_x[0] != 0);
	rev_poly(b, _solve_toeplitz_tmp[0], n);
	remainder_of_power(fast_poly_multiplication(rev_poly(_ad_x, _solve_toeplitz_tmp[1], n + 1), _solve_toeplitz_tmp[0], _solve_toeplitz_tmp[2]), n + 1);
	remainder_of_power(fast_poly_multiplication(_ad_y, _solve_toeplitz_tmp[0], _solve_toeplitz_tmp[3]), n + 1);
	return multipy_poly_by_const(remainder_of_power(sub_poly(
		fast_poly_multiplication(_ad_x, rev_poly(_solve_toeplitz_tmp[3], _solve_toeplitz_tmp[4], n), _solve_toeplitz_tmp[1]),
		fast_poly_multiplication(rev_poly(_ad_y, _solve_toeplitz_tmp[3], n + 1),
								rev_poly(_solve_toeplitz_tmp[2], _solve_toeplitz_tmp[0], n), _solve_toeplitz_tmp[2]), dst), n + 1), inverse(_ad_x[0]));
}


std::vector<size_t>& galois_field::DFT(std::vector<size_t>& src, std::vector<size_t>& dst, size_t gen_elem) {
	// call fft
	return dst;
}

std::vector<size_t>& galois_field::IDFT(std::vector<size_t>& src, std::vector<size_t>& dst, size_t gen_elem) {
	DFT(src, dst, gen_elem);
	for (size_t i = 0; i < dst.size(); ++i) {
		dst[i] = multiply(dst[i], _inverse_element);
	}
	return dst;
}

void galois_field::split_poly(std::vector<size_t> const& p, std::vector<size_t>& dst1, std::vector<size_t>& dst2, size_t m) {
	for (size_t i = 0; i < p.size(); ++i) {
		if (i < m) {
			dst1.push_back(p[i]);
		}
		else
		{
			dst2.push_back(p[i]);
		}
	}
}

std::vector<size_t>& galois_field::add_poly(std::vector<size_t>& a, std::vector<size_t>& b, std::vector<size_t>& dst, size_t m) {
	for (size_t i = 0; i < b.size() + m; ++i) {
		if (i < m && i >= a.size()) {
			a.push_back(0);
		}
		else if (i >= m) {
			if (i >= a.size()) {
				a.push_back(b[i - m]);
			}
			else {
				a.push_back(add(a[i], b[i - m]));
			}
		}
	}
	return a;
}

std::vector<size_t>& galois_field::sub_poly(std::vector<size_t>& a, std::vector<size_t>& b, std::vector<size_t>& dst) {
	size_t i;
	for (i = 0; i < a.size() && i < b.size(); ++i) {
		dst[i] = add(a[i], multiply(b[i], _p - 1));
	}
	while (i < a.size()) {
		dst[i] = a[i];
		++i;
	}
	while (i < b.size()) {
		dst[i] = multiply(b[i], _p - 1);
		++i;
	}
	return dst;
}

std::vector<size_t>& galois_field::multipy_poly_by_const(std::vector<size_t>& p, size_t a) {
	for (size_t i = 0; i < p.size(); ++i) {
		p[i] = multiply(p[i], a);
	}
	return p;
}

std::vector<size_t>& galois_field::rev_poly(std::vector<size_t>& src, std::vector<size_t>& dst, size_t n) {
	for (size_t i = 0; i < n; ++i) {
		dst[n - i - 1] = src[i];
	}
	return dst;
}

std::vector<size_t> galois_field::inv_poly(std::vector<size_t>& src, std::vector<size_t>& dst, size_t mod) {
	// make zeros in inverse tmp
	size_t g0 = inverse(src[0]);
	size_t r = std::countl_zero(mod); // mb bad
	dst[0] = 1;
	for (size_t i = 1; i < r; ++i) {
		fast_poly_multiplication(dst, dst, _inverse_temporary1);
		fast_poly_multiplication(src, _inverse_temporary1, _inverse_temporary2);
		multipy_poly_by_const(dst, 2);
		sub_poly(dst, _inverse_temporary2, _inverse_temporary1);
		using std::swap;
		swap(dst, _inverse_temporary1);
	}

	return dst;
}

std::vector<size_t>& galois_field::remainder_of_power(std::vector<size_t>& a, size_t n) {
	for (size_t i = n; i < a.size(); ++i) {
		a[i] = 0;
	}
	return a;
}


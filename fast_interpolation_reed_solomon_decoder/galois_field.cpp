#include "galois_field.h"
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <vector>
#include <bit>

#ifdef __GNUC__
#define COUNT_LEADING(v) __builtin_ctzll(v)
#else
#define COUNT_LEADING(v) _tzcnt_u64(v)
#endif

#ifdef __GNUC__
#define COUNT_TRAILING(v) __builtin_clzll(v)
#else
#define COUNT_TRAILING(v) _lzcnt_u64(v)
#endif


galois_field::galois_field(size_t p, size_t m, std::vector<size_t> const& gen) 
	: _p(p)
	, _m(m)
	, _q(std::pow(p, m))
	, _n(_q - 1)
	, _generating_poly(gen)
	, _log_table(_q) 
	, _exp_table(_q)
	, _minus_one(_p - 1)
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

size_t galois_field::poly_to_num(std::vector<size_t> const& a) const {
	size_t res = 0;
	for (auto i : a) {
		res *= _p;
		res += i;
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
	if (m < n) {
		m = _log_table[increment(a)];
		size_t log_sum = (m + n) % _n;
		return _exp_table[log_sum];
	}
	else {
		n = _log_table[increment(b)];
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

std::vector<size_t>& galois_field::fast_poly_multiplication(std::vector<size_t>& a, std::vector<size_t>& b) {
	DFT(a, _exp_table[1]);
	DFT(b, _exp_table[1]);
	for (size_t i = 0; i < a.size(); ++i) {
		a[i] = multiply(a[i], b[i]);
	}
	return IDFT(a, inverse(_exp_table[1]));
}

std::vector<size_t>& galois_field::fast_poly_division(std::vector<size_t>& a, std::vector<size_t>& b) {
	if (a.size() < b.size()) {
		return a;
	}
	//return rev_poly();
}

std::vector<std::pair<std::vector<size_t>, std::vector<size_t>>> galois_field::EMGCD(std::vector<size_t> const& a, std::vector<size_t> const& b) {
	if (b.size() < a.size() / 2) {
		return { {a, b}, {{1}, {0}}, {{0}, {1}} };
	}
	std::vector<size_t> b0, b1, c0, c1;
	size_t m = a.size() / 2 + a.size() % 2;
	std::tie(c0, b0) = split_poly(a, m);
	std::tie(c1, b1) = split_poly(b, m);

	auto u1w1v1 = EMGCD(b0, b1);
	auto u1w1v1c = u1w1v1;
	auto c0c = c0;
	auto c0cc = c0;
	
	std::vector<size_t> d,e;
	
	e = add_poly(add_poly(fast_poly_multiplication(u1w1v1[1].second, c0), fast_poly_multiplication(u1w1v1[2].second, c1), 0), u1w1v1[0].second, m);

	d = add_poly(add_poly(fast_poly_multiplication(u1w1v1[1].first, c0), fast_poly_multiplication(u1w1v1[2].first, c1), 0), u1w1v1[0].first, m);

	if (e.size() < a.size() / 2) {
		return { {d, e}, u1w1v1c[1], u1w1v1c[2] };
	}
	else {
		// q = d / e // fast division needs to be implemented
		// f = d % e
		// k = 2m - deg(e)
		// E = g0*x^k + h0
		// F = g1*x^k + h1
		// u2w2v2 = EMGCD(g0, g1)
		// return 
	}

	// assume immplemented

}

std::pair<std::vector<size_t>, std::vector<size_t>> galois_field::AD(std::vector<size_t>& a) {
	

	return {};

}

size_t galois_field::increment(size_t a) const {
	size_t r = a % _p;
	a = (a - r) + ((r + 1) % _p);
	return a;
}

std::vector<size_t>& galois_field::DFT(std::vector<size_t>& a, size_t gen_elem) {
	std::vector<size_t> a0(a.size() / 2), a1(a.size() / 2);
	for (size_t i = 0; i < a.size() / 2; ++i) {
		a0[i] = a[2 * i];
		a1[i] = a[2 * i + 1];
	}
	size_t new_gen = multiply(gen_elem, gen_elem);
	DFT(a0, new_gen);
	DFT(a1, new_gen);
	size_t base = 1;
	for (size_t i = 0; i < a.size() / 2; ++i) {
		a[i] = add(a0[i], multiply(base, a1[i]));
		a[i + a.size() / 2] = sub(a0[i], multiply(base, a1[i]));
		base = multiply(base, gen_elem);
	}
	return a;
}

std::vector<size_t>& galois_field::IDFT(std::vector<size_t>& a, size_t gen_elem) {
	DFT(a, gen_elem);
	for (size_t i = 0; i < a.size(); ++i) {
		a[i] = multiply(a[i], _inverse_element);
	}
	return a;
}

std::pair<std::vector<size_t>, std::vector<size_t>> galois_field::split_poly(std::vector<size_t> const& p, size_t m) {
	std::vector<size_t> a, b;
	for (size_t i = 0; i < p.size(); ++i) {
		if (i < m) {
			a.push_back(p[i]);
		}
		else
		{
			b.push_back(p[i]);
		}
	}
	return { std::move(a), std::move(b) };
}

std::vector<size_t>& galois_field::add_poly(std::vector<size_t>& a, std::vector<size_t>& b, size_t m) {
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



std::vector<size_t>& galois_field::rev_poly(std::vector<size_t>& a) {
	std::reverse(a.begin(), a.end());
	return a;
}

std::vector<size_t> galois_field::inv_poly(std::vector<size_t>& a, size_t mod) {
	size_t g0 = inverse(a[0]);
	size_t r = std::countl_zero(mod);



}


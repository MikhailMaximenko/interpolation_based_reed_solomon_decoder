#include "galois_field.h"
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <vector>

//#include <cstddef>


galois_field::galois_field(size_t p, size_t m, std::vector<size_t> const& gen) 
	: _p(p)
	, _m(m)
	, _q(std::pow(p, m))
	, _n(_q - 1)
	, _generating_poly(gen)
	, _log_table(_q) 
	, _exp_table(_q)
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

void galois_field::shift_poly(std::vector<size_t>& a) const {
	for (ptrdiff_t i = a.size() - 1; i >= 1; --i) {
		a[i] = a[i - 1];
	}
	a[0] = 0;
}

size_t galois_field::add(size_t, size_t) const { return 0; }
size_t galois_field::sub(size_t, size_t) const { return 0; }

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
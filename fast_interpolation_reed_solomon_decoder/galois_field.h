#pragma once
#include <cstddef>
#include<vector>


struct galois_field
{
	size_t _p;
	size_t _m;
	size_t _q;
	size_t _n;

	std::vector<size_t> _generating_poly;

	std::vector<size_t> _log_table;
	std::vector<size_t> _exp_table;

	galois_field(size_t, size_t, std::vector<size_t> const&);

	size_t multiply(size_t, size_t) const;
	size_t add(size_t, size_t) const;
	size_t sub(size_t, size_t) const;
	size_t divide(size_t, size_t) const;
	size_t inverse(size_t) const;

private:
	void init();

	std::vector<size_t> inverse_add_init(std::vector<size_t> const&) const;
	std::vector<size_t> add_init(std::vector<size_t> const&, std::vector<size_t> const&) const;
	std::vector<size_t> multiply_by_const_init(std::vector<size_t> const&, size_t) const;

	size_t poly_to_num(std::vector<size_t> const&) const;
	void shift_poly(std::vector<size_t>&) const;


};
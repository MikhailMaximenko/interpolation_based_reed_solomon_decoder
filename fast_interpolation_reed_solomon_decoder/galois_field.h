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
	size_t _minus_one;
	size_t _inverse_element; // n^(-1) (n from idft)

	// if _n is not power of 2 need to store another field over which dft/idft will be applied

	galois_field(size_t, size_t, std::vector<size_t> const&);

	size_t multiply(size_t, size_t) const;
	size_t add(size_t, size_t) const;
	size_t sub(size_t, size_t) const;
	size_t divide(size_t, size_t) const;
	size_t inverse(size_t) const;

	void fast_poly_multiplication(std::vector<size_t>&, std::vector<size_t>&);


private:
	void init();

	std::vector<size_t> inverse_add_init(std::vector<size_t> const&) const;
	std::vector<size_t> add_init(std::vector<size_t> const&, std::vector<size_t> const&) const;
	std::vector<size_t> multiply_by_const_init(std::vector<size_t> const&, size_t) const;

	size_t poly_to_num(std::vector<size_t> const&) const;
	void shift_poly(std::vector<size_t>&) const;

	size_t increment(size_t) const;

	void DFT(std::vector<size_t> &, size_t);
	void IDFT(std::vector<size_t>&, size_t);
};
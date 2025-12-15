#pragma once
#include <cstddef>
#include<vector>
#include<array>


struct galois_field
{
	// original field size
	size_t _p;

	// power of extension
	size_t _m;

	// extention field size
	size_t _q; 
	
	// multiplicative group size
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

	std::vector<size_t>& DFT(std::vector<size_t>&, std::vector<size_t>&, size_t);
	std::vector<size_t>& IDFT(std::vector<size_t>&, std::vector<size_t>&, size_t);

	std::vector<size_t>& fast_poly_multiplication(std::vector<size_t>&, std::vector<size_t>&, std::vector<size_t>&);
	std::vector<size_t>& fast_poly_division(std::vector<size_t>&, std::vector<size_t>&, std::vector<size_t>&, std::vector<size_t>&);
	std::array<std::pair<std::vector<size_t>, std::vector<size_t>>, 3>& 
		EMGCD(std::vector<size_t> const&, std::vector<size_t> const&, std::array<std::pair<std::vector<size_t>, std::vector<size_t>>, 3>&, size_t);
	void AD(std::vector<size_t>&, size_t, std::vector<size_t>&, std::vector<size_t>&);

	std::vector<size_t>& SOLVE_TOEPITZ(std::vector<size_t>&, std::vector<size_t>&, size_t, std::vector<size_t>&);

	size_t degree(std::vector<size_t> const&);
	std::vector<size_t>& rev_poly(std::vector<size_t>&, std::vector<size_t>&, size_t);
	std::vector<size_t>& add_poly(std::vector<size_t>&, std::vector<size_t>&, std::vector<size_t>&, size_t);
	std::vector<size_t>& sub_poly(std::vector<size_t>&, std::vector<size_t>&, std::vector<size_t>&);

private:
	void init();

	std::vector<size_t> inverse_add_init(std::vector<size_t> const&) const;
	std::vector<size_t> add_init(std::vector<size_t> const&, std::vector<size_t> const&) const;
	std::vector<size_t> multiply_by_const_init(std::vector<size_t> const&, size_t) const;

	size_t poly_to_num(std::vector<size_t> const&) const;
	std::vector<size_t>& shift_poly(std::vector<size_t>&) const;


	void split_poly(std::vector<size_t> const&, std::vector<size_t>&, std::vector<size_t>&, size_t);

	std::vector<size_t>& multipy_poly_by_const(std::vector<size_t>&, size_t);
	

	std::vector<size_t>& remainder_of_power(std::vector<size_t>&, size_t);
	std::vector<size_t> inv_poly(std::vector<size_t>&, std::vector<size_t>&, size_t);

private:
	std::vector<size_t> _a_tmp;
	std::vector<size_t> _b_tmp;
	std::vector<size_t> _inverse_temporary1;
	std::vector<size_t> _inverse_temporary2;
	std::vector<std::array<std::pair<std::vector<size_t>, std::vector<size_t>>, 3>> _emgcd_tmp_result;
	std::vector<std::array<std::vector<size_t>, 8>> _emgcd_tmp_polynomials;

	std::array<std::array<std::pair<std::vector<size_t>, std::vector<size_t>>, 3>, 2> _ad_tmp_emgcd_results;
	std::array<std::vector<size_t>, 3> _ad_tmp_polinomyals;

	std::array<std::vector<size_t>, 12> _solve_toeplitz_tmp;

	std::vector<size_t> _ad_x;
	std::vector<size_t> _ad_y;
	std::vector<size_t> _b_rev;

};
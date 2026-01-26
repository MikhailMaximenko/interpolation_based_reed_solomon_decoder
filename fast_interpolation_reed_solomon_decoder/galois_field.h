#pragma once
#include <cstddef>
#include<vector>
#include<array>


struct galois_field
{
	// original field size
	//unsigned _p;

	// power of extension
	unsigned _m;

	// extention field size
	unsigned _q; 
	
	// multiplicative group size
	unsigned _n; 

	unsigned _gen_poly;
	unsigned _poly_size;
	std::vector<unsigned> _exp_table;
	std::vector<unsigned> _log_table;

	unsigned _inverse_element; // n^(-1) (n from idft)

	// if _n is not power of 2 need to store another field over which dft/idft will be applied

	galois_field(unsigned, unsigned, unsigned);

	unsigned multiply(unsigned, unsigned) const;
	unsigned add(unsigned, unsigned) const;
	unsigned divide(unsigned, unsigned) const;
	unsigned inverse(unsigned) const;

	std::vector<unsigned>& DFT(std::vector<unsigned>&, std::vector<unsigned>&); // cyclotomic
	std::vector<unsigned>& IDFT(std::vector<unsigned>&, std::vector<unsigned>&); // cyclotomic

	std::vector<unsigned>& DFT(std::vector<unsigned>&, std::vector<unsigned>&, unsigned, unsigned); // binary architecture a with variable size
	std::vector<unsigned>& IDFT(std::vector<unsigned>&, std::vector<unsigned>&, unsigned, unsigned); // binary architecture a with variable size

	std::vector<unsigned>& fast_poly_multiplication(std::vector<unsigned>&, std::vector<unsigned>&, std::vector<unsigned>&);
	std::vector<unsigned>& fast_poly_division(std::vector<unsigned>&, std::vector<unsigned>&, std::vector<unsigned>&, std::vector<unsigned>&);
	std::array<std::pair<std::vector<unsigned>, std::vector<unsigned>>, 3>& 
		EMGCD(std::vector<unsigned> const&, std::vector<unsigned> const&, std::array<std::pair<std::vector<unsigned>, std::vector<unsigned>>, 3>&, unsigned);
	void AD(std::vector<unsigned>&, unsigned, std::vector<unsigned>&, std::vector<unsigned>&);

	std::vector<unsigned>& SOLVE_TOEPITZ(std::vector<unsigned>&, std::vector<unsigned>&, unsigned, std::vector<unsigned>&);

	unsigned degree(std::vector<unsigned> const&);
	std::vector<unsigned>& rev_poly(std::vector<unsigned>&, std::vector<unsigned>&, unsigned);
	std::vector<unsigned>& add_poly(std::vector<unsigned>&, std::vector<unsigned>&, std::vector<unsigned>&, unsigned);
	std::vector<unsigned>& sub_poly(std::vector<unsigned>&, std::vector<unsigned>&, std::vector<unsigned>&);

private:
	void init();

	std::vector<unsigned>& DFTimpl(std::vector<unsigned>&, std::vector<unsigned>&, unsigned, unsigned, unsigned); // binary architecture a with variable size
	std::vector<unsigned>& IDFTimpl(std::vector<unsigned>&, std::vector<unsigned>&, unsigned, unsigned, unsigned); // binary architecture a with variable size

	//std::vector<unsigned> inverse_add_init(std::vector<unsigned> const&) const;
	//std::vector<unsigned> add_init(std::vector<unsigned> const&, std::vector<unsigned> const&) const;
	//std::vector<unsigned> multiply_by_const_init(std::vector<unsigned> const&, unsigned) const;

	unsigned poly_to_num(std::vector<unsigned> const&) const;
	std::vector<unsigned>& shift_poly(std::vector<unsigned>&) const;


	void split_poly(std::vector<unsigned> const&, std::vector<unsigned>&, std::vector<unsigned>&, unsigned);

	std::vector<unsigned>& multipy_poly_by_const(std::vector<unsigned>&, unsigned);
	

	std::vector<unsigned>& remainder_of_power(std::vector<unsigned>&, unsigned);
	std::vector<unsigned> inv_poly(std::vector<unsigned>&, std::vector<unsigned>&, unsigned);

public:
	std::vector<unsigned> _a_tmp;
	std::vector<unsigned> _b_tmp;
	std::vector<unsigned> _inverse_temporary1;
	std::vector<unsigned> _inverse_temporary2;
	std::vector<std::array<std::pair<std::vector<unsigned>, std::vector<unsigned>>, 3>> _emgcd_tmp_result;
	std::vector<std::array<std::vector<unsigned>, 8>> _emgcd_tmp_polynomials;

	std::array<std::array<std::pair<std::vector<unsigned>, std::vector<unsigned>>, 3>, 2> _ad_tmp_emgcd_results;
	std::array<std::vector<unsigned>, 3> _ad_tmp_polinomyals;

	std::array<std::vector<unsigned>, 12> _solve_toeplitz_tmp;

	std::vector<std::vector<std::vector<unsigned>>> _dft_tmp;

	std::vector<std::vector<unsigned>> _s;

	std::vector<unsigned> _ad_x;
	std::vector<unsigned> _ad_y;
	std::vector<unsigned> _b_rev;

};
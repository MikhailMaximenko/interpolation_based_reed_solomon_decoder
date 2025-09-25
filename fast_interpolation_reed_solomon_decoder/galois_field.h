#pragma once
#include<vector>


struct galois_field
{
	size_t _q;
	size_t _m;

	std::vector<size_t> _generating_poly;

	std::vector<size_t> _log_table;
	std::vector<size_t> _exp_table;

	galois_field(size_t, size_t, std::vector<size_t> const&);

	size_t multiply(size_t, size_t) const;
	size_t add(size_t, size_t) const;
	size_t sub(size_t, size_t) const;
	size_t divide(size_t, size_t) const;




};
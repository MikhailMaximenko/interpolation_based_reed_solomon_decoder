#include "galois_field.h"

//#include <cstddef>


galois_field::galois_field(size_t q, size_t m, std::vector<size_t> const& gen) 
	: _q(q)
	, _m(m)
	, _generating_poly(gen)
{


	
}

size_t galois_field::multiply(size_t, size_t) const { return 0; }
size_t galois_field::add(size_t, size_t) const { return 0; }
size_t galois_field::sub(size_t, size_t) const { return 0; }
size_t galois_field::divide(size_t, size_t) const { return 0; }
// fast_interpolation_reed_solomon_decoder.cpp : Defines the entry point for the application.
//

#include "fast_interpolation_reed_solomon_decoder.h"


#ifdef __GNUC__
#define COUNT_LEADING(v) __builtin_ctzll(v)
#else
#include<intrin.h>
#define COUNT_LEADING(v) __tzcnt(v)
#endif

#ifdef __GNUC__
#define COUNT_TRAILING(v) __builtin_clzll(v)
#else
//#include<intrin.h>
#define COUNT_TRAILING(v) __lzcnt(v)
#endif

int main()
{
	uint64_t a = 8;
	std::cout << COUNT_LEADING(a) << "\n";
	return 0;
}

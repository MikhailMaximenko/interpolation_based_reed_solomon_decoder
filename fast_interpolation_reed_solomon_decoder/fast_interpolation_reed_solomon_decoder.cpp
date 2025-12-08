// fast_interpolation_reed_solomon_decoder.cpp : Defines the entry point for the application.
//

#include "fast_interpolation_reed_solomon_decoder.h"
#include "galois_field.h"
#include <bit>

int main()
{
	galois_field gf(11, 1, { 9, 1 });
	for (auto i : gf._log_table) {
		std::cout << i << " ";
	}
	std::cout << "\n";
	for (auto i : gf._exp_table) {
		std::cout << i << " ";
	}
	std::cout << "\n";

	std::cout << gf.divide(8, 10) << "\n";


	return 0;
}

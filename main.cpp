#include "fft.h"
#include "berlecamp_massey_decoder.h"
#include "fast_interpolation_reed_solomon_decoder.h"
#include <cassert>
#include <bit>
#include <fstream>
#include <random>

std::random_device rd;
std::mt19937 gen(rd());


std::vector<unsigned> generate_message(unsigned n, unsigned k) {
	std::vector<unsigned> msg(n);
	std::uniform_int_distribution<unsigned> distr(0, n);
	std::generate_n(msg.begin(), k, [&]() { return distr(gen); });
	return msg;
}

std::vector<unsigned> generate_errors(unsigned n, unsigned t) {
	std::vector<unsigned> errors(n);
	std::uniform_int_distribution<unsigned> distr(1, n);
	std::generate_n(errors.begin(), t, [&]() { return distr(gen); });
	std::shuffle(errors.begin(), errors.end(), gen);
	return errors;
}


template<typename DecType>
void test_decoder(DecType& decoder, unsigned n, unsigned k, unsigned t, unsigned iters, std::ofstream& metrics_dump) {
	//encoding::bch_decoder decoder(gf, n, k, 2, n - k + 1, 1);
	//InterpolationBasedFastRSDecoder decoder(gf, n, k);
	//for (size_t t = (n - k ) / 2 - 10; t <= (n - k)/2; ++t) {
	std::cout << "testing t: " << t << "\n";
	for (size_t _ = 0; _ < iters; ++_) {
		auto msg = generate_message(n, k);
		std::vector<unsigned> encoded = decoder.encode(msg);
		auto errors = generate_errors(n, t);
		std::vector<unsigned> msg_with_errors(n);
		decoder._gf.add_poly(encoded, errors, msg_with_errors, 0);
		decoder._gf.reset_counters();
		decoder.decode(msg_with_errors);

		std::cout << "additions: " << decoder._gf._additions << " multiplications: " << decoder._gf._multiplications << "\n";
		metrics_dump << n << " " << k << " " << t << " " << decoder._gf._additions << " " << decoder._gf._multiplications << "\n";
		decoder._gf.reset_counters();
		for (size_t i = 0; i < n; ++i) {
			if (encoded[i] != msg_with_errors[i]) {
				std::cout << "decoding error occured with t:" << t << "\n";
				decoder._gf.print_poly(msg);
				decoder._gf.print_poly(encoded);
				decoder._gf.print_poly(errors);
				decoder._gf.print_poly(msg_with_errors);
				return;
			}
		}
		//}
		//break;
	}
	std::cout << "tests passed\n";
}

void benchmark_multiplication(std::string out_prefix, unsigned method) {
	std::vector<unsigned> field_sizes = { 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16 };
	std::vector<unsigned> field_generators = { 0xb, 0x13, 0x25, 0x43, 0x83, 0x11d, 0x211, 0x409, 0x805, 0x1053, 0x201b, 0x4143, 0x8003, 0x1100b };
	if (method == 0) {
		unsigned initial = 3;
		for (size_t i = 0; i < 5; ++i) {
			std::string o_name = out_prefix + std::to_string(initial) + std::string(".txt");
			std::ofstream o(o_name);
			for (size_t j = 0; j < 8; ++j) {
				unsigned q = 1 << field_sizes[j];
				galois_field gf(field_sizes[j], field_generators[j], field_sizes[j]);
				gf._naive_mult_schonhage = initial;
				gf.reset_counters();
				std::vector<unsigned> a(2 * q, 1), b(2 * q, 1), c(2 * q);
				unsigned len = q - 1;
				unsigned n = 1;
				while (len >= 1) {
					n *= 3;
					len /= 3;
				}
				gf.SCHONHAGE_STRASSEN_FFT(a, b, c, n, 0);
				o << q << " " << gf._additions << " " << gf._multiplications << "\n";
			}
			initial *= 3;
		}
	}
	else if (method == 1) {
		unsigned initial = 2;
		for (size_t i = 0; i < 7; ++i) {
			std::string o_name = out_prefix + std::to_string(initial) + std::string(".txt");
			std::ofstream o(o_name);
			for (size_t j = 0; j < 8; ++j) {
				unsigned q = 1 << field_sizes[j];
				galois_field gf(field_sizes[j], field_generators[j], field_sizes[j]);
				gf._naive_mult_karatsuba= initial;
				gf.reset_counters();
				std::vector<unsigned> a(q, 1), b(q, 1), c(2 * q);
				
				gf.karatsuba_multiplication(a, b, c, q, 0);
				o << q << " " << gf._additions << " " << gf._multiplications << "\n";
			}
			initial *= 2;
		}
	}
	else if (method == 2) {
		unsigned initial = 9;
		for (size_t i = 0; i < 7; ++i) {
			std::string o_name = out_prefix + std::to_string(initial) + std::string(".txt");
			std::ofstream o(o_name);
			for (size_t j = 0; j < 9; ++j) {
				unsigned q = 1 << field_sizes[j];
				galois_field gf(field_sizes[j], field_generators[j], field_sizes[j]);
				gf._karatsuba_mult_schonhage = initial;
				gf.reset_counters();
				std::vector<unsigned> a(2 * q, 1), b(2 * q, 1), c(2 * q);
				unsigned len = q - 1;
				unsigned n = 1;
				while (len >= 1) {
					n *= 3;
					len /= 3;
				}
				gf.SCHONHAGE_STRASSEN_FFT(a, b, c, n, 0);
				o << q << " " << gf._additions << " " << gf._multiplications << "\n";
			}
			initial *= 3;
		}
	}
	else if (method == 3) {
		unsigned initial = 1;
		std::string o_name = out_prefix + std::to_string(initial) + std::string(".txt");
		std::ofstream o(o_name);
		for (size_t j = 0; j < 11; ++j) {
			unsigned q = 1 << field_sizes[j];
			galois_field gf(field_sizes[j], field_generators[j], field_sizes[j]);
			gf._karatsuba_mult_gao_mateer = initial;
			gf.reset_counters();
			std::vector<unsigned> a(2 * q, 1), b(2 * q, 1), c(2 * q), d(2 * q), e(2 * q), f(2 * q);
			unsigned m = field_sizes[j];
			gf.gao_mateer_fft(a, c, m);
			gf.gao_mateer_fft(b, d, m);
			for (size_t k = 0; k < q; ++k) {
				e[k] = gf.multiply(c[k], d[k]);
			}
			gf.gao_mateer_ifft(e, f, m);
			o << q / 2 << " " << gf._additions << " " << gf._multiplications << "\n";
		}
		//initial *= 3;
	}


}

// 1: output
// 2: code rate
// 3: error rate
// 4: decoder type: 0 - interpolation decoder, 1 - berlecamp-massey decoder
int main(int argc, char* argv[])
{
	//double k_rate, err_rate;
	//int type;
	//std::ofstream out("interpolation_results_0_5rate_max_errs_long.txt");
	////out << "hello\n";
	//if (false) {
	//	if (argc != 4) {
	//		std::cout << "expected 3 args: output file, code rate (from 0 to 1), error rate (from 0 to 1)\n";
	//		return -1;
	//	}
	//	try {
	//		k_rate = std::stod(argv[2]);
	//		err_rate = std::stod(argv[3]);
	//		type = std::stoi(argv[4]);
	//	}
	//	catch (std::invalid_argument const& e) {
	//		std::cout << "could not parse numeric args: " << e.what();
	//		return -1;
	//	}
	//}
	//else {
	//	k_rate = 0.5;
	//	err_rate = 0.9;
	//	type = 0;
	//}
	//if (out.bad()) {
	//	std::cout << "could not open input file\n";
	//	return -1;
	//}

	//out << "n, k, err num, additions, multiplications\n";
	////////galois_field gf2(7, 0x83, 7);

	////////test_decoder(gf2, 127, 11, 15);

	//////galois_field gf(3, 0xb, 3);
	//std::vector<unsigned> field_sizes = { 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16};
	//std::vector<unsigned> field_generators = { 0xb, 0x13, 0x25, 0x43, 0x83, 0x11d, 0x211, 0x409, 0x805, 0x1053, 0x201b, 0x4143, 0x8003, 0x1100b };
	////// 2 2
	//for (size_t i = 0; i < field_sizes.size(); ++i) {
	//	galois_field gf(field_sizes[i], field_generators[i], field_sizes[i]);
	//	unsigned n = (1 << field_sizes[i]) - 1;
	//	unsigned k = n * k_rate;
	//	unsigned t = ((n - k) / 2) * err_rate;
	//	if (type == 0) {
	//		InterpolationBasedFastRSDecoder decoder(std::move(gf), n, k);
	//		test_decoder(decoder, n, k, t, 1, out);
	//	}
	//	else if (type == 1) {
	//		encoding::bch_decoder decoder(std::move(gf), n, k, 2, n - k + 1, 1);
	//		std::cout << "here\n";
	//		test_decoder(decoder, n, k, t, 1, out);
	//		//test_decoder(decoder, 511, 256, t, 10, out);
	//	}
	//	else {
	//		std::cout << "unknown decoder type\n";
	//		return -1;
	//	}
	//}

	benchmark_multiplication("gao_mateer", 3);
	return 0;

	
}
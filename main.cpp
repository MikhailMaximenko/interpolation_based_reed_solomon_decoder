#include "fft.h"
#include "berlecamp_massey_decoder.h"
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

// 1: output
// 2: code rate
// 3: error rate
// 4: decoder type: 0 - interpolation decoder, 1 - berlecamp-massey decoder
int main(int argc, char* argv[])
{
	//double k_rate, err_rate;
	//int type;
	//std::ofstream out("berlecamp_decoder_results_0_75rate_half_errs.txt");
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
	//	k_rate = 0.75;
	//	err_rate = 0.5;
	//	type = 1;
	//}
	//if (out.bad()) {
	//	std::cout << "could not open input file\n";
	//	return -1;
	//}

	//out << "n, k, err num, additions, multiplications\n";
	//////galois_field gf2(7, 0x83, 7);

	//////test_decoder(gf2, 127, 11, 15);

	//////galois_field gf(3, 0xb, 3);
	//std::vector<unsigned> field_sizes = { 3, 6, 7, 9, 10, 11 };
	//std::vector<unsigned> field_generators = { 0xb, 0x43, 0x83, 0x211, 0x409, 0x805 };
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
	//return 0;

	
}
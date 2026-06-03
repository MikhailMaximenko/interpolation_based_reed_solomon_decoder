#include "galois_field.h"
#include <gtest/gtest.h>

TEST(taylor_expansion, one_quater_poly) {
	galois_field gf(3, 0xb, 3);

	std::vector<unsigned > a = { 1, 1 };
	std::vector<unsigned> b(2), c(2);
	gf.taylor_expansion(a, b, 2, 2);
	std::vector<unsigned> check = { 1, 1 };
	EXPECT_EQ(b, check);
	gf.itaylor_expansion(b, c, 2, 2);
	EXPECT_EQ(a, c);
}

TEST(taylor_expansion, random_quater_poly) {
	galois_field gf(3, 0xb, 3);

	std::vector<unsigned > a = { 3, 4 };
	std::vector<unsigned> b(2), c(2);
	gf.print_poly(gf.precomputed_basises_delta[1]);
	gf.taylor_expansion(a, b, 2, 2);
	std::vector<unsigned> check = { 3, 4 };
	EXPECT_EQ(b, check);
	gf.itaylor_expansion(b, c, 2, 2);
	EXPECT_EQ(a, c);
}

TEST(taylor_expansion, one_half_poly) {
	galois_field gf(3, 0xb, 3);

	std::vector<unsigned > a = { 1, 1, 1, 1, 0, 0, 0, 0 };
	std::vector<unsigned> b(8), c(8);
	gf.taylor_expansion(a, b, 4, 2);
	std::vector<unsigned> check = { 1, 1, 0, 1, 0, 0, 0, 0 };
	EXPECT_EQ(b, check);
	gf.itaylor_expansion(b, c, 4, 2);
	EXPECT_EQ(a, c);
}

TEST(taylor_expansion, random_half_poly) {
	galois_field gf(3, 0xb, 3);

	std::vector<unsigned > a = { 3, 3, 2, 7, 0, 0, 0, 0 };
	std::vector<unsigned> b(8), c(8);
	gf.taylor_expansion(a, b, 4, 2);
	std::vector<unsigned> check = { 3, 6, 5, 7, 0, 0, 0, 0 };
	EXPECT_EQ(b, check);
	gf.itaylor_expansion(b, c, 4, 2);
	EXPECT_EQ(a, c);
}

TEST(taylor_expansion, full_one_poly) {
	galois_field gf(3, 0xb, 3);

	std::vector<unsigned > a = { 1, 1, 1, 1, 1, 1, 1, 1 };
	std::vector<unsigned> b(8), c(8);
	gf.taylor_expansion(a, b, 8, 2);
	std::vector<unsigned> check = { 1, 1, 0, 1, 0, 0, 0, 1 };
	EXPECT_EQ(b, check);
	gf.itaylor_expansion(b, c, 8, 2);
	gf.print_poly(c);
	EXPECT_EQ(a, c);
}

TEST(taylor_expansion, full_random_poly) {
	galois_field gf(3, 0xb, 3);

	std::vector<unsigned > a = { 2, 4, 4, 5, 1, 3, 3, 6 };
	std::vector<unsigned> b(8), c(8);
	gf.taylor_expansion(a, b, 8, 2);
	std::vector<unsigned> check = { 2, 2, 6, 0, 2, 0, 5, 6 };
	EXPECT_EQ(b, check);
	gf.itaylor_expansion(b, c, 8, 2);
	gf.print_poly(c);
	EXPECT_EQ(a, c);
}

TEST(init_galois_field, gao_basises_check) {
	galois_field gf(3, 0xb, 3);

	std::vector<unsigned> delta3 = { 1, 2, 4 };
	std::vector<unsigned> gamma2 = { 7, 5};
	std::vector<unsigned> delta2 = { 4, 2 };
	std::vector<unsigned> gamma1 = { 2 };
	std::vector<unsigned> delta1 = { 6 };
	EXPECT_EQ(delta3, gf.precomputed_basises_delta[3]);
	EXPECT_EQ(delta2, gf.precomputed_basises_delta[2]);
	EXPECT_EQ(delta1, gf.precomputed_basises_delta[1]);
	EXPECT_EQ(gamma2, gf.precomputed_basises_gamma[2]);
	EXPECT_EQ(gamma1, gf.precomputed_basises_gamma[1]);
}

TEST(init_galois_field, gao_gamma_spaces) {
	galois_field gf(3, 0xb, 3);

	std::vector<unsigned> gamma_space2 = {0, 7, 5, 2 };
	std::vector<unsigned> gamma_space1 = {0, 2 };
	EXPECT_EQ(gamma_space2, gf.precomputed_space_gamma[2]);
	EXPECT_EQ(gamma_space1, gf.precomputed_space_gamma[1]);
}

TEST(gao_mateer_fft, one_quarter_poly) {
	galois_field gf(3, 0xb, 3);
	std::vector<unsigned> a = { 1, 1 };
	std::vector<unsigned> b(2), c(2);
	gf.gao_mateer_fft(a, b, 1);
	std::vector<unsigned> check = {
		gf.substitute_poly(a, 0),
		gf.substitute_poly(a, 6),
	};
	EXPECT_EQ(b, check);

	gf.gao_mateer_ifft(b, c, 1);
	EXPECT_EQ(a, c);
}

TEST(gao_mateer_fft, random_quarter_poly) {
	galois_field gf(3, 0xb, 3);
	std::vector<unsigned> a = { 6, 2 };
	std::vector<unsigned> b(2), c(2);
	gf.gao_mateer_fft(a, b, 1);
	std::vector<unsigned> check = {
		gf.substitute_poly(a, 0),
		gf.substitute_poly(a, 6),
	};
	EXPECT_EQ(b, check);

	gf.gao_mateer_ifft(b, c, 1);
	EXPECT_EQ(a, c);
}

TEST(gao_mateer_fft, one_half_poly) {
	galois_field gf(3, 0xb, 3);
	std::vector<unsigned> a = {1, 1, 1, 1};
	std::vector<unsigned> b(4), c(4);
	gf.gao_mateer_fft(a, b, 2);
	std::vector<unsigned> check = {
		gf.substitute_poly(a, 0),
		gf.substitute_poly(a, 4),
		gf.substitute_poly(a, 2),
		gf.substitute_poly(a, 6)
	};
	EXPECT_EQ(b, check);

	gf.gao_mateer_ifft(b, c, 2);
	EXPECT_EQ(a, c);
}

TEST(gao_mateer_fft, random_half_poly) {
	galois_field gf(3, 0xb, 3);
	std::vector<unsigned> a = { 2, 7, 6, 2};
	std::vector<unsigned> b(4), c(4);
	gf.gao_mateer_fft(a, b, 2);
	std::vector<unsigned> check = {
		gf.substitute_poly(a, 0),
		gf.substitute_poly(a, 4),
		gf.substitute_poly(a, 2),
		gf.substitute_poly(a, 6)
	};
	EXPECT_EQ(b, check);
	gf.gao_mateer_ifft(b, c, 2);
	EXPECT_EQ(a, c);
}


TEST(gao_mateer_fft, one_full_poly) {
	galois_field gf(3, 0xb, 3);
	std::vector<unsigned> a = { 1, 1, 1, 1, 1, 1, 1, 1 };
	std::vector<unsigned> b(8), c(8);

	gf.gao_mateer_fft(a, b, 3);
	std::vector<unsigned> check = {
		gf.substitute_poly(a, 0),
		gf.substitute_poly(a, 1),
		gf.substitute_poly(a, 2),
		gf.substitute_poly(a, 3),
		gf.substitute_poly(a, 4),
		gf.substitute_poly(a, 5),
		gf.substitute_poly(a, 6),
		gf.substitute_poly(a, 7),
	};
	EXPECT_EQ(b, check);
	gf.gao_mateer_ifft(b, c, 3);
	EXPECT_EQ(a, c);
}


TEST(gao_mateer_fft, random_full_poly) {
	galois_field gf(3, 0xb, 3);
	std::vector<unsigned> a = { 7, 3, 2, 2, 1, 6, 7, 3 };
	std::vector<unsigned> b(8), c(8);

	gf.gao_mateer_fft(a, b, 3);
	std::vector<unsigned> check = {
		gf.substitute_poly(a, 0),
		gf.substitute_poly(a, 1),
		gf.substitute_poly(a, 2),
		gf.substitute_poly(a, 3),
		gf.substitute_poly(a, 4),
		gf.substitute_poly(a, 5),
		gf.substitute_poly(a, 6),
		gf.substitute_poly(a, 7),
	};
	EXPECT_EQ(b, check);
	gf.gao_mateer_ifft(b, c, 3);
	EXPECT_EQ(a, c);
}




TEST(gao_mateer_fft, random_full_poly_large_field) {
	galois_field gf(6, 0x43, 6);
	std::vector<unsigned> a = { 15, 24, 13, 61, 46, 21, 61, 11, 32, 33, 23, 15, 41, 31, 32, 3  };
	std::vector<unsigned> b(16), c(16);

	gf.gao_mateer_fft(a, b, 4);
	gf.gao_mateer_ifft(b, c, 4);
	EXPECT_EQ(a, c);
}

TEST(gao_mateer_fft, multiply_short_poly) {
	galois_field gf(3, 0xb, 3);
	std::vector<unsigned> 
		a = { 2, 3, 0, 0 }, 
		b = { 2, 3, 0, 0 }, 
		check = { 4, 0, 5, 0};
	std::vector<unsigned> c(4), d(4), e(4), f(4);
	gf.gao_mateer_fft(a, c, 2);
	gf.gao_mateer_fft(b, d, 2);
	for (size_t i = 0; i < 4; ++i) {
		e[i] = gf.multiply(c[i], d[i]);
	}
	gf.gao_mateer_ifft(e, f, 2);
	EXPECT_EQ(f, check);
}

TEST(gao_mateer_fft, multiply_short_long_poly) {
	galois_field gf(3, 0xb, 3);
	std::vector<unsigned>
		a = { 2, 3, 0, 0, 0, 0, 0, 0 },
		b = { 2, 3, 6, 0, 0, 0, 0, 0 },
		check = { 4, 0, 2, 1, 0, 0, 0, 0 };
	std::vector<unsigned> c(8), d(8), e(8), f(8);
	gf.gao_mateer_fft(a, c, 3);
	gf.gao_mateer_fft(b, d, 3);
	for (size_t i = 0; i < 8; ++i) {
		e[i] = gf.multiply(c[i], d[i]);
	}
	gf.gao_mateer_ifft(e, f, 3);
	EXPECT_EQ(f, check);
}


TEST(gao_mateer_fft, operations_check) {
	std::vector<unsigned> field_sizes = { 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14};
	std::vector<unsigned> field_generators = { 0xb, 0x13, 0x25, 0x43, 0x83, 0x11d, 0x211, 0x409, 0x805, 0x1053, 0x201b, 0x4143};


	for (size_t i = 0; i < field_sizes.size(); ++i) {
		unsigned m = field_sizes[i];
		galois_field gf(m, field_generators[i], m);

		std::vector<unsigned> a(1 << m), b(1 << m);
		gf.gao_mateer_fft(a, b, m);
		gf.gao_mateer_fft(a, b, m);
		gf.gao_mateer_ifft(a, b, m);
		gf.reset_counters();
	}

}
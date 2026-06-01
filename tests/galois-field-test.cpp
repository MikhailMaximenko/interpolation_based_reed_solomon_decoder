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
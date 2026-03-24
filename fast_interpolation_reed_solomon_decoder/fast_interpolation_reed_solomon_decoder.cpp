// fast_interpolation_reed_solomon_decoder.cpp : Defines the entry point for the application.
//

#include "fast_interpolation_reed_solomon_decoder.h"
#include "galois_field.h"
#include "fft.h"
#include <bit>

InterpolationBasedFastRSDecoder::InterpolationBasedFastRSDecoder(galois_field const& gf, unsigned n, unsigned k)
	: _gf(gf)
	, _n(n)
	, _k(k)
	, _t((n - k) / 2)
{
	for (auto& v : _tmp) {
		v.resize(_n + 1);
	}
}


std::vector<unsigned> InterpolationBasedFastRSDecoder::encode(std::vector<unsigned> & m) {
	std::vector<unsigned> res(_n);
	_gf.DFT(m, res);
	return res;
}


void InterpolationBasedFastRSDecoder::decode(std::vector<unsigned>& cw) {
	//_gf.print_poly(cw);
	//_gf.print_poly(_tmp[0]);
	//_t = 0;
	if (_t != 0) {
		for (auto& tmp : _tmp) {
			std::fill(tmp.begin(), tmp.end(), 0);
		}
		_gf.IDFT(cw, _tmp[0]);
		std::copy(_tmp[0].begin() + _n - 2 * _t + 1, _tmp[0].begin() + _n, _tmp[1].begin());
		std::copy(_tmp[0].begin() + _n - 2 * _t, _tmp[0].begin() + _n - _t, _tmp[2].begin());
		//std::copy(_tmp[1].begin(), _tmp[1].begin() + 2 * _t - 1, _tmp[0].rbegin());
		//std::copy(_tmp[2].begin(), _tmp[2].begin() + _t, _tmp[0].begin() + _n - 2 * _t);
		//std::cout << _t << " initial system:\n";
		for (auto& v : _tmp[1]) {
			//std::cout << _gf._log_table[v] << " ";
		}
		//std::cout << "\n";
		for (auto& v : _tmp[2]) {
			//std::cout << _gf._log_table[v] << " ";
		}
		//std::cout << "\n";
		//_gf.print_poly(_tmp[1]);
		//_gf.print_poly(_tmp[2]);
		//_tmp[1].insert(_tmp[1].begin(), _tmp[0].rbegin(), _tmp[0].rend() - (_n - 2 * _t));
		//_tmp[2].insert(_tmp[2].begin(), _tmp[0].rbegin() - _t, _tmp[0].rend() - (_n - 2 * _t - 1)); // n - 2 * t + 1 ..n - t
		_gf.SOLVE_TOEPITZ(_tmp[1], _tmp[2], _t - 1, _tmp[3]);
		//std::cout << "solutions:\n";
		//_gf.print_poly(_tmp[3]);
		for (auto& v : _tmp[3]) {
			//std::cout << _gf._log_table[v] << " ";
		}
		//std::cout << "\n";

		for (ptrdiff_t i = _k - 1; i >= 0; --i) {
			_tmp[0][i] = 0;
			for (size_t j = 1; j <= _t; ++j) {
				_tmp[0][i] = _gf.add(_tmp[0][i], _gf.multiply(_tmp[0][i + j], _tmp[3][_t - j]));
			}
			//std::cout << _tmp[0][i] << " ";
		}
		//std::cout << "\n";
		//_gf.print_poly(_tmp[0]);

		using std::swap;
		swap(cw, _tmp[1]);
		//_gf.print_poly(_tmp[0]);
		_gf.DFT(_tmp[0], _tmp[2]);
		//_gf.IDFT(_tmp[2], _tmp[0]);
		//_gf.print_poly(_tmp[0]);
		_gf.sub_poly(_tmp[1], _tmp[2], cw);
		//_gf.print_poly(cw);
	}



}

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



void test_decoder(galois_field & gf, unsigned n, unsigned k, unsigned iters) {
	InterpolationBasedFastRSDecoder decoder(gf, n, k);
	for (size_t t = 0; t <= (n - k) / 2; ++t) {
		std::cout << "testing t: " << t << "\n";
		decoder._t = t;
		for (size_t _ = 0; _ < iters; ++_) {
			//std::cout << "here\n";
			auto msg = generate_message(n, k);
			auto encoded = decoder.encode(msg);
			auto errors = generate_errors(n, t);
			std::vector<unsigned> msg_with_errors(n + 1);
			gf.add_poly(encoded, errors, msg_with_errors, 0);
			//gf.print_poly(msg_with_errors);
			decoder.decode(msg_with_errors);
			for (size_t i = 0; i < n; ++i) {
				if (encoded[i] != msg_with_errors[i]) {
					std::cout << "decoding error occured with t:"<< t << "\n";
					gf.print_poly(msg);
					gf.print_poly(encoded);
					gf.print_poly(errors);
					gf.print_poly(msg_with_errors);
					return;
				}
			}
		}
	}
	std::cout << "tests passed\n";


}

int main()
{
	galois_field gf2(7, 0x83, 7);
	/*std::vector <unsigned> a1{0, 2, 4, 0, 0, 0, 0, 0};
	std::vector <unsigned> b1{ 4, 0, 0, 0, 0, 0, 0, 0 };
	std::vector<unsigned> ds{ 0, 0, 0, 0, 0, 0, 0, 0 };
	gf2.SOLVE_TOEPITZ(a1, b1, 1, ds);
	gf2.print_poly(ds);*/
	/*
	std::cout << "here\n";*/
	for (size_t i = 0; i < (1 << 6); ++i) {
		std::cout << gf2._exp_table[i] << " ";
	}
	std::cout << "\n";
	for (size_t i = 0; i < (1 << 6); ++i) {
		std::cout << gf2._log_table[i] << " ";
	}

	std::cout << "\n";
	/*
	std::vector<unsigned> m1 = { 3, 2, 0, 0, 0, 0, 0, 0 }, m2 = { 1, 1, 0, 0, 0, 0, 0, 0 }, m3 = { 0,0,0,0,0,0,0,0}, m4 = { 0,0,0,0,0,0,0,0}, m5 = { 0,0,0,0,0,0,0,0 }, m6 = { 0,0,0,0,0,0,0,0 };
	gf2.DFT(m1, m3, 8, 0);
	gf2.DFT(m2, m4, 8, 0);
	for (size_t i = 0; i < 8; ++i) {
		m5[i] = gf2.multiply(m3[i], m4[i]);
	}

	gf2.print_poly(m3);
	gf2.print_poly(m4);
	gf2.print_poly(m5);
	gf2.IDFT(m5, m6, 8, 0);*/
	//gf2.multipy_poly_by_const(m6, 3);
	/*for (auto v : m6) {
		std::cout << v << " ";
	}
	std::cout << "\n";*/



	//gf2.fast_poly_multiplication(m1, m2, m3);
	//gf2.print_poly(m3);

	//gf2.inv_poly(m1, m3, 7);
	//gf2.print_poly(m3);

	/*std::vector<unsigned> a{ 1, 0, 0, 1, 0, 0, 0, 0 };
	std::vector<unsigned> b{ 0, 1, 0, 0, 0, 0, 0, 0 };*/
	//std::vector<unsigned> e{ 3, 6, 1, 0, 0, 0, 0, 0 };
	//std::vector<unsigned> f{ 0, 0, 0, 1, 0, 0, 0, 0 };
	////
	////std::vector<unsigned> c(8), d(8);
	////gf2.SOLVE_TOEPITZ(a, b, 2, c);
	/////*gf2.AD(a, 2, c, d);
	////gf2.print_poly(d);*/
	////gf2.print_poly(c);

	////std::cout << "!\n";
	/*std::array<std::pair<std::vector<unsigned>, std::vector<unsigned> >, 3> dst =
	{
	  {
		{ {0,0,0,0,0,0,0,0}, {0,0,0,0,0,0,0,0} },
		{ {0,0,0,0,0,0,0,0}, {0,0,0,0,0,0,0,0} },
		{ {0,0,0,0,0,0,0,0}, {0,0,0,0,0,0,0,0} }
	  }
	};
	std::cout << "calling emgcd\n";
	gf2.EMGCD(f, e, dst, 0);
	for (auto& v : dst) {
		gf2.print_poly(v.first);
		gf2.print_poly(v.second);
	}*/

	//call_fft(gf2._b_tmp, gf2._a_tmp);
	//for (auto i : gf2._a_tmp) {
		//std::cout << i << " ";
	//}
	//std::cout << "\n";
	//std::cout << "res:\n";
	//for (auto i : c) {
	//	std::cout << i << " ";
	//}
	//std::cout << "\n";


	//std::cout << "snd\n";
	//std::vector<unsigned> a1 = { 2, 5, 0, 0, 0, 0, 0 }, b1 = {2, 5, 0, 0, 0, 0, 0}, f, g;
	//f.resize(7); g.resize(7);
	//gf2.DFT(a1, f);
	//std::cout << "fst\n";
	//gf2.DFT(b1, g);

	//std::cout << "completed!\n";

	//for (size_t i = 0; i < a1.size(); ++i) {
	//	a1[i] = gf2.multiply(f[i], g[i]);
	//	std::cout << a1[i] << " ";
	//}
	//std::cout << "\n";

	//gf2.IDFT(a1, b1);
	//for (auto i : b1) {
	//	std::cout << i << " ";
	//}
	//std::cout << "\n";
	//a1 = { 2, 3, 0, 0, 0, 0, 0, 0 };
	//gf2.DFT(a1, f, 4, 0);
	//std::cout << "straight\n";
	//for (auto i : f) {
	//	std::cout << i << " ";
	//}
	//std::cout << "\n";
	//std::cout << "inverse\n";
	//for (auto i : gf2.IDFT(f, g, 4, 0)) {
	//	std::cout << i << " ";
	//}
	//std::cout << "\n";

	//std::cout << "id expected for:\n";
	//a1 = { 2, 1, 3, 5, 2, 7, 3, 1 };
	//b1 = { 0, 0, 0, 0, 0, 0, 0, 0 };
	//for (auto i : a1) {
	//	std::cout << i << " ";
	//}
	//std::cout << "\ndft:\n";
	//gf2.DFT(a1, b1, 8, 0);
	//for (auto i : b1) {
	//	std::cout << i << " ";
	//}
	//for (auto& i : a1) {
	//	i = 0;
	//}
	//std::cout << "\n";
	//gf2.IDFT(b1, a1, 8, 0);
	//std::cout << "got:\n";
	//for (auto i : a1) {
	//	std::cout << i << " ";
	//}
	//std::cout << "\n";

	//std::vector<unsigned> aa{ 2, 5, 0, 0, 0, 0, 0, 0 };
	//std::vector<unsigned> bb{ 2, 5, 0, 0, 0, 0, 0, 0 };
	//std::vector<unsigned> cc{ 4, 0, 7, 0, 0, 0, 0, 0 };

	//std::vector<unsigned> r1(8), r2(8), r3(8);

	//gf2.DFT(aa, r1, 4, 0);
	//std::cout << "dft1:\n";
	//for (auto i : r1) {
	//	std::cout << i << " ";
	//}
	//std::cout << "\n";
	//gf2.DFT(bb, r2, 4, 0);
	//std::cout << "dft2:\n";
	//for (auto i : r2) {
	//	std::cout << i << " ";
	//}
	//std::cout << "\n";
	//for (size_t i = 0; i < 8; ++i) {
	//	r3[i] = gf2.multiply(r1[i], r2[i]);
	//	r1[i] = 0;
	//}
	//for (auto i : r3) {
	//	std::cout << i << " ";
	//}
	//std::cout << "\n";
	//gf2.IDFT(r3, r1, 4, 0);
	//for (auto kk : r1) {
	//	std::cout << kk << " ";
	//}
	//std::cout << "\n";
	//gf2.DFT(cc, r2, 4, 0);
	//std::cout << "expected result dft:\n";
	//for (auto i : r2) {
	//	std::cout << i << " ";
	//}
	//std::cout << "\n";
	/*gf2.IDFT(r1, r2, 4, 0);
	std::cout << "idft1:\n";
	for (auto i : r2) {
		std::cout << i << " ";
	}
	std::cout << "\n";*/
	/*std::vector<unsigned> a_inv(7), kk(7);
	gf2.inv_poly(a, a_inv, 7);

	gf2.fast_poly_multiplication(a, a_inv, kk);
	for (auto v : a_inv) {
		std::cout << v << " ";
	}
	std::cout << "\n";
	for (auto v : kk) {
		std::cout << v << " ";
	}
	std::cout << "\n";
	*/


	//std::vector<unsigned> a_inv(7), kk(7);

	//gf2.inv_poly(a, a_inv, 7);

	//gf2.fast_poly_multiplication(a, a_inv, kk);

	//for (auto i : a_inv) {
	//	std::cout << i << " ";
	//}
	//std::cout << "\n";
	//for (auto i : kk) {
	//	std::cout << i << " ";
	//}
	//std::cout << "\n";


	//std::vector<unsigned> aaa = { 2, 3, 0, 0, 0, 0, 0};
	//std::vector<unsigned> ccc(7), bbb = { 1, 0, 0, 0, 0, 0, 0};
	//gf2.fast_poly_multiplication(aaa, bbb, ccc);
	//gf2.print_poly(ccc);

	//std::vector<unsigned> m = {1,2,3,0,0,0,0,0};
	//InterpolationBasedFastRSDecoder decoder(gf2, 63, 5);
	auto a = gf2._ad_tmp_polinomyals[0];
	auto e = gf2._ad_tmp_polinomyals[0];
	auto b = a;
	auto c = a;
	auto d = c;
	auto f = c;
	auto dst = gf2._emgcd_tmp_result[0]; //42 25 56 42 8 48 10
	a[0] = 42;
	a[1] = 25;
	a[2] = 56;
	a[3] = 42;
	a[4] = 8;
	a[5] = 48;
	a[6] = 10;
	gf2.shift_poly(a);
	a[8] = 1;
	a[0] = 1;
	b[0] = 19;
	b[1] = 20;
	b[2] = 28;
	e[9] = 1;

	//59 53 52 39 57 45 53 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
//46 59 53 52 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0

	//gf2.SOLVE_TOEPITZ(a, b, 2, c);
	//gf2.print_poly(c);
	//std::cout << "-----\n";
	////gf2.AD(a, 3, c, d);
	////gf2.print_poly(c);
	////gf2.print_poly(d);
	//std::cout << "-----\n";
	/*gf2.EMGCD(e, a, dst, 0);
	for (auto& v : dst) {
		gf2.print_poly(v.first);
		gf2.print_poly(v.second);
	}*/
	//std::cout << "-----\n";
	//gf2.rev_poly(a,f,8);
	//////gf2.print_poly(f);
	//for (auto& v : dst) {
	//	std::fill(v.first.begin(), v.first.end(), 0);
	//	std::fill(v.second.begin(), v.second.end(), 0);
	//}
	////gf2.print_poly(e);
	////gf2.print_poly(f);
	//gf2.EMGCD(e, f, dst, 0);
	//for (auto& v : dst) {
	//	gf2.print_poly(v.first);
	//	gf2.print_poly(v.second);
	//}

	test_decoder(gf2, 127, 11, 15);
	return 0;
}
//53 45 29 53 3 10 13 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
//23 53 45 29 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
//42 25 56 42 8 48 10 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
//41 42 25 56 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
//solutions :
//43 44 33 20 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
//51 37 62 14 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
//decoding error occured with t : 4
//4 27 13 6 53 8 37 21 19 14 20 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
//16 58 6 9 30 58 9 9 46 59 1 1 9 23 11 2 12 23 29 56 35 48 43 45 1 22 22 53 53 19 41 9 28 4 44 54 9 2 55 27 44 6 9 31 11 18 61 52 5 26 25 30 43 19 9 48 56 51 63 10 56 20 28
//0 40 0 0 0 0 0 62 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 57 0 0 0 0 0 0 0 36 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
//41 27 12 20 55 55 17 53 16 12 27 41 13 5 35 14 15 14 48 61 21 56 2 47 58 13 29 55 3 12 2 48 18 29 48 13 40 35 29 35 27 28 48 21 39 36 4 61 6 18 46 9 52 60 63 20 26 58 47 34 57 24 63 0
// fast_interpolation_reed_solomon_decoder.cpp : Defines the entry point for the application.
//

#include "fast_interpolation_reed_solomon_decoder.h"
#include "galois_field.h"
#include "fft.h"
#include <bit>



void InterpolationBasedFastRSDecoder::decode(std::vector<unsigned>& cw) {
	_gf.IDFT(cw, _tmp[0]);
	_tmp[1].insert(_tmp[1].begin(), _tmp[0].rbegin(), _tmp[0].rend() + (_n - 2 * _t));
	_tmp[2].insert(_tmp[2].begin(), _tmp[0].rbegin() - _t, _tmp[0].rend() + (_n - 2 * _t - 1)); // n - 2 * t + 1 ..n - t
	_gf.SOLVE_TOEPITZ(_tmp[1], _tmp[2], _t, _tmp[3]);
	for (ptrdiff_t i = _k - 1; i >= 0; --i) {
		_tmp[0][i] = 0;
		for (size_t j = 0; j < _t; ++j) {
			_tmp[0][i] = _gf.add(_tmp[0][i], _gf.multiply(_tmp[0][i + j], _tmp[3][j]));
		}
	}

	using std::swap;
	swap(cw, _tmp[1]);
	_gf.sub_poly(_tmp[1], _gf.DFT(_tmp[0], _tmp[2]), cw);


}

int main()
{
	galois_field gf2(3, 0xb, 3);
	std::cout << "here\n";
	for (size_t i = 0; i < (1 << 3); ++i) {
		std::cout << gf2._exp_table[i] << " ";
	}
	std::cout << "\n";
	for (size_t i = 0; i < (1 << 3); ++i) {
		std::cout << gf2._log_table[i] << " ";
	}

	std::cout << "\n";

	std::vector<unsigned> a{ 0, 2, 1, 7, 0, 0, 0 };
	std::vector<unsigned> b{ 0, 2, 3, 0, 0, 0, 0 };
	
	std::vector<unsigned> c(7), d(7);
	std::array<std::pair<std::vector<unsigned>, std::vector<unsigned> >, 3> dst = 
	{ 
      {
	    { {0,0,0,0,0,0,0}, {0,0,0,0,0,0,0} }, 
		{ {0,0,0,0,0,0,0}, {0,0,0,0,0,0,0} }, 
		{ {0,0,0,0,0,0,0}, {0,0,0,0,0,0,0} }
	  }
	};
	std::cout << "calling emgcd\n";
	gf2.EMGCD(a, b, dst, 0);
	for (auto v : dst) {
		gf2.print_poly(v.first);
	}
	for (auto v : dst) {
		gf2.print_poly(v.second);
	}
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


	return 0;
}

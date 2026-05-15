#include "berlecamp_massey_decoder.h"
#include "linalg.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <vector>

namespace encoding {

 
    bch_decoder::bch_decoder(galois_field fld, size_t n, size_t k, size_t root, size_t delta, size_t b)
        : _n(n)
        , _k(k)
        , _gf(std::move(fld))
        , _root(root)
        , _delta(delta)
        , _b(b)
        , _decoding(_n * fld._m)
        , _helper_poly(_delta - 1, 0)
        , _tmp_poly(_delta - 1, 0)
        , _syndromes(_delta - 1, 0)
        , _locators_poly(_n, 0)
        , _translation(_n)
        , _reordering(n)
        , _locators_values(n)
        , _locators_inversed_roots(n)
        , _locators_roots(n)
        , _locators_derivate(n)
        , _error_values(n)
        , _derivate_values(n)
        , _forney_values(n)
        , _substitution_poly(n)
        , _forney(n)
    {
    }

    std::vector<unsigned> bch_decoder::encode(std::vector<unsigned>& cw) {
        std::vector<unsigned> res(_n);
        _gf.DFT(cw, res);
        return res;
    }

    size_t bch_decoder::berlecamp_massey(std::vector<unsigned>& syndroms) {
        _helper_poly[0] = 1;
        _locators_poly[0] = 1;
        size_t l, m;
        m = 0;
        l = 0;
        for (size_t r = 1; r <= _delta - 1; ++r) {
            size_t Delta = 0;
            for (size_t j = 0; j <= l; ++j) {
                Delta = _gf.add(Delta, _gf.multiply(_locators_poly[j], syndroms[r - 1 - j]));
            }

            if (Delta != 0) {
                for (size_t i = 0; i < r - m; ++i) {
                    _tmp_poly[i] = _locators_poly[i];
                }
                for (size_t i = r - m; i < _delta - 1; ++i) {
                    _tmp_poly[i] = _gf.add(_locators_poly[i], _gf.multiply(Delta, _helper_poly[i - (r - m)]));
                }

                if (2 * l <= r - 1) {
                    for (size_t i = 0; i < _delta - 1; ++i) {
                        _helper_poly[i] = _gf.multiply(_gf.inverse(Delta), _locators_poly[i]);
                    }
                    l = r - l;
                    m = r;
                }
                std::copy_n(_tmp_poly.begin(), _delta - 1, _locators_poly.begin());
            }

        }
        return l;
    }



    std::vector<unsigned>& bch_decoder::pgz(std::vector<unsigned>& cw) {
        size_t cur_beta = _gf._exp_table[_b];
        size_t start_beta = cur_beta;
        bool flag = true;

        std::fill_n(_locators_poly.begin(), _n, 0);
        std::fill_n(_translation.begin(), _n, 0);
        std::fill_n(_reordering.begin(), _n, 0);
        std::fill_n(_locators_values.begin(), _n, 0);
        std::fill_n(_locators_inversed_roots.begin(), _n, 0);
        std::fill_n(_locators_roots.begin(), _n, 0);
        std::fill_n(_locators_derivate.begin(), _n, 0);
        std::fill_n(_error_values.begin(), _n, 0);
        std::fill_n(_derivate_values.begin(), _n, 0);
        std::fill_n(_forney_values.begin(), _n, 0);
        std::fill_n(_substitution_poly.begin(), _n, 0);
        std::fill_n(_forney.begin(), _n, 0);

        //_gf.translate_bit_vector(_decoding, _translation);

        //std::cout << _decoding.to_string() << "\n";

        //_gf.print_poly(_translation);

        // count syndromes 
        _gf.DFT(cw, _reordering);
        //std::cout << _gf._log_table[_root] << " " << start_beta << "\n";
        for (size_t i = 0; i < _delta - 1; ++i) {
            //std::cout << _gf._log_table[cur_beta] << " ";
            size_t cur = 1;

            _syndromes[i] = _reordering[_gf._log_table[cur_beta]];
            if (_syndromes[i]) {
                flag = false;
            }
            cur_beta = _gf.multiply(cur_beta, _root);
        }
        //std::cout << "\n";
        //_gf.print_poly(_reordering);
        // std::cout << "syndromes: ";
        //_gf.print_poly(_syndromes);
        // for (auto i : _syndromes) {
            // std::cout << i << " ";
        // }
        // std::cout << "\n";
        if (flag) {
            return cw;
        }

        auto sz = berlecamp_massey(_syndromes);
        for (ptrdiff_t b = _locators_poly.size() - 1; b > 0; --b) {
            if (_locators_poly[b] != 0) {
                if (b == sz) {
                    break;
                }
                else {
                    return cw;
                }
            }
        }

        //std::cout << "locators: ";
        _gf.DFT(_locators_poly, _locators_values);
        //_gf.print_poly(_locators_poly);
        //_gf.print_poly(_locators_values);
        size_t roots_cnt = 0;
        for (size_t i = 0; i < _gf._n; ++i) {
            if (_locators_values[i] == 0) {
                size_t pos = _gf.inverse(_gf._exp_table[i]);
                // std::cout << pos << " ";
                _locators_roots[roots_cnt] = _gf._exp_table[i];
                _locators_inversed_roots[roots_cnt] = pos;
                ++roots_cnt;
            }
        }
        //std::cout << "locators roots: ";
        //_gf.print_poly(_locators_roots);
        //for (size_t v : _locators_roots) {
        //    std::cout << _gf._log_table[v] << " ";
        //}

        //std::cout << std::endl;
        //_gf.print_poly(_locators_inversed_roots);
        //for (size_t v : _locators_inversed_roots) {
        //    std::cout << _gf._log_table[v] << " ";
        //}

         std::cout << std::endl;
         size_t t = (_delta - 1) / 2;
        _gf.remainder_of_power(_gf.fast_poly_multiplication(_locators_poly, _syndromes, _forney), 2 * t);
        //_gf.formal_derivate(_locators_poly, _locators_derivate);
        _gf.DFT(_forney, _forney_values);
        _gf.DFT(_locators_derivate, _derivate_values);
        for (size_t i = 0; i < roots_cnt; ++i) {
            unsigned root_pos = _gf._log_table[_locators_roots[i]] % _gf._n;
            unsigned val = _gf.multiply(_locators_roots[i], _forney_values[root_pos]);
            //std::cout << val << " " << _locators_roots[i] << " " << _forney_values[root_pos] << "\n";
            for (size_t j = 0; j < roots_cnt; ++j) {
                if (i != j) {
                    val = _gf.multiply(val, _gf.inverse(_gf.add(1, _gf.multiply(_locators_roots[i], _locators_inversed_roots[j]))));
                }
            }
            //std::cout << val << "\n";
            cw[_gf._log_table[_locators_inversed_roots[i]] % _gf._n] = _gf.add(cw[_gf._log_table[_locators_inversed_roots[i]] % _gf._n], val);
        }
        //_gf.print_poly(_translation);
        _gf.translate_to_bit_vector(_translation, _decoding);

        for (size_t i = 0; i < _delta - 1; ++i) {
            _tmp_poly[i] = 0;
            _locators_poly[i] = 0;
            _helper_poly[i] = 0;
        }
        return cw;
    }

    std::vector<unsigned>& bch_decoder::decode(std::vector<unsigned>& cw) {
        return pgz(cw);
    }

    linalg::bit_vector bch_decoder::decode() {
        _gf.translate_bit_vector(_decoding, _translation);
        // std::cout << _decoding.to_string() << "\n";
        pgz(_translation);

        // std::cout << _decoding.to_string() << "\n";
        // linalg::bit_vector decoded;
        // for (auto const& vect : _gen) {
        //     if (vect.leading() == _decoding.leading()) {
        //         decoded.push_back(true);
        //         _decoding += vect;
        //     } else {
        //         decoded.push_back(false);
        //     }
        // }
        // std::cout << _decoding.to_string() << "\n";
        return linalg::bit_vector(_decoding);
    }

    linalg::bit_vector bch_decoder::decode(std::vector<double> const& signals) {
        std::vector<double> reliabilities;
        linalg::bit_vector hard_decisions;
        for (size_t i = 0; i < signals.size(); ++i) {
            if (signals[i] >= 0.0) {
                _decoding.set(i, true);
            }
            else {
                _decoding.set(i, false);
            }
        }
        return decode();
    }

}
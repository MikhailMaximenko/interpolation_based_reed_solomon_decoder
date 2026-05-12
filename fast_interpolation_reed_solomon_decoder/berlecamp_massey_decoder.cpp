#include "berlecamp_massey_decoder.h"
#include "linalg.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <vector>

namespace encoding {

 
    bch_decoder::bch_decoder(galois_field const& fld, size_t n, size_t k, size_t root, size_t delta, size_t b)
        : _n(n)
        , _k(k)
        , _fld(fld)
        , _root(_fld._exp_table[root])
        , _delta(delta)
        , _b(b)
        , _decoding(_n)
        , _locators_poly(_n, 0)
        , _helper_poly(_delta - 1, 0)
        , _tmp_poly(_delta - 1, 0)
        , _syndromes(_delta - 1, 0)
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


    size_t bch_decoder::berlecamp_massey(std::vector<unsigned>& syndroms) {
        _helper_poly[0] = 1;
        _locators_poly[0] = 1;
        size_t l, m;
        m = 0;
        l = 0;
        for (size_t r = 1; r <= _delta - 1; ++r) {
            size_t Delta = 0;
            for (size_t j = 0; j <= l; ++j) {
                Delta = _fld.add(Delta, _fld.multiply(_locators_poly[j], syndroms[r - 1 - j]));
            }

            if (Delta != 0) {
                for (size_t i = 0; i < r - m; ++i) {
                    _tmp_poly[i] = _locators_poly[i];
                }
                for (size_t i = r - m; i < _locators_poly.size(); ++i) {
                    _tmp_poly[i] = _fld.add(_locators_poly[i], _fld.multiply(Delta, _helper_poly[i - (r - m)]));
                }

                if (2 * l <= r - 1) {
                    for (size_t i = 0; i < _locators_poly.size(); ++i) {
                        _helper_poly[i] = _fld.multiply(_fld.inverse(Delta), _locators_poly[i]);
                    }
                    l = r - l;
                    m = r;
                }
                std::copy_n(_tmp_poly.begin(), _delta - 1, _locators_poly.begin());
            }

        }
        return l;
    }



    void bch_decoder::pgz() {
        size_t cur_beta = _root;
        size_t start_beta = cur_beta;
        bool flag = true;
        _fld.translate_bit_vector(_decoding, _translation);

        // count syndromes 
        _fld.DFT(_translation, _reordering);
        for (size_t i = 0; i < _delta - 1; ++i) {
            size_t cur = 1;

            _syndromes[i] = _reordering[cur_beta];
            if (_syndromes[i]) {
                flag = false;
            }
            cur_beta = _fld.multiply(cur_beta, start_beta);
        }
        // std::cout << "syndromes: ";
        // for (auto i : _syndromes) {
            // std::cout << i << " ";
        // }
        // std::cout << "\n";
        if (flag) {
            return;
        }

        auto sz = berlecamp_massey(_syndromes);
        for (ptrdiff_t b = _locators_poly.size() - 1; b > 0; --b) {
            if (_locators_poly[b] != 0) {
                if (b == sz) {
                    break;
                }
                else {
                    return;
                }
            }
        }

        std::cout << "locators: ";
        _fld.DFT(_locators_poly, _locators_values);
        size_t roots_cnt = 0;
        for (size_t i = 0; i < _fld._n; ++i) {
            if (_locators_values[i] == 0) {
                size_t pos = _fld.inverse(i) % (_fld._n);
                // std::cout << pos << " ";
                _locators_roots[roots_cnt] = i;
                _locators_inversed_roots[roots_cnt] = pos;
                ++roots_cnt;
            }
        }
        // std::cout << std::endl;
        size_t t = _delta / 2;
        _fld.remainder_of_power(_fld.fast_poly_multiplication(_locators_poly, _syndromes, _forney), 2 * t + 1);
        _fld.formal_derivate(_locators_poly, _locators_derivate);
        _fld.DFT(_forney, _forney_values);
        _fld.DFT(_locators_derivate, _derivate_values);
        for (size_t i = 0; i < roots_cnt; ++i) {
            unsigned root_pos = _fld._log_table[_locators_roots[i]] % _fld._n;
            unsigned val = _fld.multiply(_fld.multiply(_locators_inversed_roots[i], _forney_values[root_pos]), _fld.inverse(_derivate_values[root_pos]));
            _translation[_fld._log_table[_locators_inversed_roots[i]] % _fld._n] = _fld.add(_translation[_fld._log_table[_locators_inversed_roots[i]] % _fld._n], val);
        }

        for (size_t i = 0; i < _delta - 1; ++i) {
            _tmp_poly[i] = 0;
            _locators_poly[i] = 0;
            _helper_poly[i] = 0;
        }
    }

    linalg::bit_vector bch_decoder::decode() {
        // std::cout << _decoding.to_string() << "\n";
        pgz();

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
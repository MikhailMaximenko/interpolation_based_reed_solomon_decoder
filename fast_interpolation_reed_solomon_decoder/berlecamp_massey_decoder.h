#pragma once

#include "linalg.h"
#include "galois_field.h"
#include <array>
#include <cmath>
#include <cstddef>
#include <iterator>
#include <limits>
#include <memory>
#include <ostream>
#include <set>
#include <stdexcept>
#include <unordered_map>
#include <utility>
#include <vector>

namespace encoding {


    struct bch_decoder {
        galois_field _gf;
        size_t _n;
        size_t _k;
        size_t _root;
        size_t _delta;
        size_t _b;
        linalg::bit_vector _decoding;
        std::vector<unsigned> _translation;
        std::vector<unsigned> _reordering;
        std::vector<unsigned> _locators_poly;
        std::vector<unsigned> _locators_values;
        std::vector<unsigned> _locators_inversed_roots;
        std::vector<unsigned> _locators_roots;
        std::vector<unsigned> _locators_derivate;
        std::vector<unsigned> _error_values;
        std::vector<unsigned> _derivate_values;
        std::vector<unsigned> _forney_values;
        std::vector<unsigned> _helper_poly;
        std::vector<unsigned> _substitution_poly;
        std::vector<unsigned> _forney;
        std::vector<unsigned> _tmp_poly;
        std::vector<unsigned> _syndromes;


        bch_decoder(galois_field const&, size_t, size_t, size_t , size_t, size_t);

        size_t berlecamp_massey(std::vector<unsigned>&);

        void pgz();

        linalg::bit_vector decode();

        linalg::bit_vector decode(std::vector<double> const& signals);

    };

}
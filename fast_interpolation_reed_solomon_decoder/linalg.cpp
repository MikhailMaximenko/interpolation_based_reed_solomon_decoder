#include "linalg.h"
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <iostream>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace linalg {

#ifdef __GNUC__
#define COUNT_LEADING(v) __builtin_ctzll(v);
#else
#define COUNT_LEADING(v) _tzcnt_u64(v);
#endif

#ifdef __GNUC__
#define COUNT_TRAILING(v) __builtin_clzll(v);
#else
#define COUNT_TRAILING(v) _lzcnt_u64(v);
#endif

    bit_vector::bit_vector(size_t v, size_t sz) : _storage(1, v), _sz(sz) {}
    bit_vector::bit_vector() : _storage(), _sz() {}
    void bit_vector::push_back(bool b) {
        if (_sz % UNDERLIING_TYPE_SIZE == 0) {
            _storage.push_back(0);
        }
        if (b) {
            _storage[_sz / UNDERLIING_TYPE_SIZE] = _storage[_sz / UNDERLIING_TYPE_SIZE] | (1ull << (_sz % UNDERLIING_TYPE_SIZE));
        }
        ++_sz;
    }
    bit_vector::bit_vector(size_t sz, bool v) : _storage(sz / UNDERLIING_TYPE_SIZE + (sz % UNDERLIING_TYPE_SIZE ? 1 : 0), v ? -1ull : 0), _sz(sz) {
        if (v && (sz % UNDERLIING_TYPE_SIZE != 0)) {
            _storage.back() &= ((1ull << (sz % UNDERLIING_TYPE_SIZE)) - 1);
        }
    }

    bit_vector::bit_vector(std::vector<bool> const& o) : _sz(0) {
        for (size_t i = 0; i < o.size(); ++i) {
            push_back(o[i]);
        }
    }

    void bit_vector::resize(size_t sz) {
        _storage.resize(sz / UNDERLIING_TYPE_SIZE + (sz % UNDERLIING_TYPE_SIZE ? 1 : 0));
        _sz = sz;
    }

    size_t bit_vector::size() const noexcept {
        return _sz;
    }
    bool bit_vector::empty() const noexcept {
        return _sz == 0;
    }

    bool bit_vector::operator[](size_t ind) const {
        // std::cout << "op[]: " << _sz << " " << ind << "\n"; 
        return (_storage[ind / UNDERLIING_TYPE_SIZE] >> (ind % UNDERLIING_TYPE_SIZE)) & 1ull;
    }
    void bit_vector::set(size_t ind, bool v) {
        if (v) {
            _storage[ind / UNDERLIING_TYPE_SIZE] = _storage[ind / UNDERLIING_TYPE_SIZE] | (1ull << (ind % UNDERLIING_TYPE_SIZE));
        }
        else {
            _storage[ind / UNDERLIING_TYPE_SIZE] = _storage[ind / UNDERLIING_TYPE_SIZE] & ~(1ull << (ind % UNDERLIING_TYPE_SIZE));
        }
    }

    void bit_vector::cp(bit_vector const& o) {
        if (o._storage.size() > _storage.size()) {
            _storage.resize(o._storage.size());
        }
        for (size_t i = 0; i < o._storage.size(); ++i) {
            _storage[i] = o._storage[i];
        }
        _sz = o._sz;
    }

    bit_vector& bit_vector::operator+=(bit_vector const& o) {
        if (size() < o.size()) { // assume that small ctor is filled with zeros
            throw std::domain_error("cannot add vectors of smaller and bigger size");
        }
        for (size_t i = 0; i < o._storage.size(); ++i) {
            _storage[i] ^= o._storage[i];
        }
        return *this;
    }
    bit_vector bit_vector::operator+(bit_vector const& o) const {
        bit_vector tmp(*this);
        tmp += o;
        return tmp;
    }

    bit_vector& bit_vector::add_with_offset(bit_vector const& o, size_t offset) {
        if (size() < o.size() + offset) { // assume that small ctor is filled with zeros
            throw std::domain_error("cannot add vectors of smaller and bigger size");
        }
        size_t blocks = offset / UNDERLIING_TYPE_SIZE;
        size_t inner = offset % UNDERLIING_TYPE_SIZE;
        uint64_t prev = 0;
        for (size_t i = 0; i < o._storage.size(); ++i) {
            _storage[i + blocks] ^= (o._storage[i] << inner) | (prev & ((1ull << inner) - 1));
            prev = o._storage[i] >> (UNDERLIING_TYPE_SIZE - inner);
        }
        return *this;
    }


    bit_vector bit_vector::operator-() const {
        bit_vector res(size());
        for (size_t i = 0; i < _storage.size(); ++i) {
            res._storage[i] = ~_storage[i];
        }
        if (_sz % UNDERLIING_TYPE_SIZE) {
            res._storage.back() = (~_storage.back()) & ((1ull << (_sz % UNDERLIING_TYPE_SIZE)) - 1);
        }
        return res;
    }

    bit_vector& bit_vector::multiply(bit_vector const& o) {
        if (size() != o.size()) {
            throw std::domain_error("cannot multiply vectors of different size");
        }
        for (size_t i = 0; i < _storage.size(); ++i) {
            _storage[i] &= o._storage[i];
        }
        return *this;
    }

    bit_vector bit_vector::operator*(bit_matrix const& mt) const {
        if (_sz != mt.size()) {
            throw std::domain_error("operator*: unmatched sizes");
        }
        bit_vector res(mt[0].size());
        for (size_t i = 0; i < _sz; ++i) {
            if ((*this)[i]) {
                res += mt[i];
            }
        }
        return res;
    }


    bit_vector& bit_vector::operator++() noexcept {
        for (size_t i = 0; i < _storage.size(); ++i) {
            ++_storage[i];
            if (_storage[i] != 0) {
                return *this;
            }
            if (i == _storage.size() - 1 && _storage[i] == (1ull << (_sz % UNDERLIING_TYPE_SIZE))) {
                _storage[i] = 0;
                return *this;
            }
        }
        return *this;
    }

    bool bit_vector::any(bit_vector const& mask) const {
        if (mask.size() != size()) {
            throw std::domain_error("uncompatible sizes");
        }
        for (size_t i = 0; i < mask._storage.size(); ++i) {
            if ((_storage[i] & mask._storage[i]) != 0) {
                return true;
            }
        }
        return false;
    }

    bool bit_vector::all(bit_vector const& mask) const {
        if (mask.size() != size()) {
            throw std::domain_error("uncompatible sizes");
        }
        for (size_t i = 0; i < mask._storage.size(); ++i) {
            if ((_storage[i] & (~mask._storage[i])) == 0) {
                return true;
            }
        }
        return false;
    }

    size_t bit_vector::leading() const noexcept {
        size_t res = 0;
        for (size_t i = 0; i < _storage.size(); ++i) {
            if (_storage[i] == 0) {
                res += 64;
            }
            else {
                res += COUNT_LEADING(_storage[i]);
            }
            if (res != UNDERLIING_TYPE_SIZE) {
                return std::min(res, _sz);
            }
        }
        return _sz;
    }
    size_t bit_vector::trailing() const noexcept {
        for (size_t i = _sz - 1; i > 0; --i) {
            if ((*this)[i]) {
                return i;
            }
        }
        return 0;
    }

    bool bit_vector::all_zeros(size_t from, size_t to) const {
        bit_vector mask(size(), true);
        mask += bit_vector(to, true);
        if (from) mask += bit_vector(from - 1, true);
        return all(mask);
    }

    bit_vector bit_vector::puncture(size_t from, size_t to) const {
        bit_vector res(to - from);
        uint64_t shift = from % UNDERLIING_TYPE_SIZE;
        for (size_t i = 0; i < (to - from) / UNDERLIING_TYPE_SIZE; ++i) {
            res._storage[i] = (_storage[i + from / UNDERLIING_TYPE_SIZE] >> shift) | (shift ? (_storage[i + from / UNDERLIING_TYPE_SIZE + 1] << (UNDERLIING_TYPE_SIZE - shift)) : 0);
        }
        if ((to - from) % UNDERLIING_TYPE_SIZE) {
            if ((to - 1) % UNDERLIING_TYPE_SIZE < from % UNDERLIING_TYPE_SIZE) {
                res._storage.back() = (_storage[to / UNDERLIING_TYPE_SIZE - 1] >> shift) & ((1ull << (UNDERLIING_TYPE_SIZE - from % UNDERLIING_TYPE_SIZE)) - 1) |
                    ((_storage[to / UNDERLIING_TYPE_SIZE] & (((1ull << (to % UNDERLIING_TYPE_SIZE)) - 1))) << (UNDERLIING_TYPE_SIZE - from % UNDERLIING_TYPE_SIZE));
            }
            else {
                res._storage.back() = (_storage[(to - 1) / UNDERLIING_TYPE_SIZE] >> shift) & ((1ull << (to - from) % UNDERLIING_TYPE_SIZE) - 1);
            }
        }
        return res;
    }

    bit_vector bit_vector::concat(bit_vector const& o) const { // need to be fast
        bit_vector res(*this);
        res.concat(o);
        return res;
    }

    bit_vector& bit_vector::concat(bit_vector const& o) {
        if ((_storage.size() > (_sz + o._sz - 1) / UNDERLIING_TYPE_SIZE) && _storage.size() == 1) {
            _storage[0] |= o._storage[0] << _sz;
            _sz += o._sz;
            return *this;
        }
        size_t old_size = _storage.size();
        size_t old_sz = _sz;
        uint64_t shift = size() % UNDERLIING_TYPE_SIZE;
        resize(size() + o.size());
        if (shift == 0) {
            for (size_t i = old_size; i < _storage.size(); ++i) {
                _storage[i] = o._storage[i - _storage.size()];
            }
        }
        else {
            uint64_t prev = _storage.back();
            for (size_t i = old_size - 1; i < _storage.size(); ++i) {
                _storage[i] = prev;
                if (i + 1 - old_size < o._storage.size()) {
                    _storage[i] |= (o._storage[i + 1 - old_size] << shift);
                    prev = (o._storage[i + 1 - old_size] >> (UNDERLIING_TYPE_SIZE - shift));
                }
            }
        }
        return *this;
    }

    std::string bit_vector::to_string() const {
        std::string res;
        for (size_t i = 0; i < _sz; ++i) {
            res += (*this)[i] ? '1' : '0';
        }
        return res;
    }


    bool bit_vector::operator==(bit_vector const& o) const noexcept {
        return (_sz == o._sz) && (_storage == o._storage);
    }
    bool bit_vector::operator!=(bit_vector const& o) const noexcept {
        return (_sz != o._sz) || (_storage != o._storage);
    }
    bool bit_vector::operator<(bit_vector const& o) const noexcept {
        // if (size() < o.size()) {
        //     return true;
        // }
        for (ptrdiff_t i = _storage.size() - 1; i >= 0; --i) {
            if (_storage[i] < o._storage[i]) {
                return true;
            }
        }
        return false;
    }

    uint64_t bit_vector::to_bit_mask() const noexcept {
        if (_storage.empty()) {
            return 0;
        }
        return _storage[0];
    }

    std::vector<size_t> bit_matrix::get_g_s(ptrdiff_t h) const {
        if (empty()) {
            return {};
        }
        std::vector<size_t> res;
        bit_vector mask1(at(0).size(), false);
        mask1 += bit_vector(h, true);
        bit_vector mask2(at(0).size(), true);
        mask2 += bit_vector(h, true);

        for (size_t i = 0; i < size(); ++i) {
            bit_vector const& v = at(i);
            if (v.any(mask2) && v.any(mask1)) {
                res.push_back(i);
            }
        }
        return res;
    }

    std::vector<size_t> bit_matrix::get_g_f(ptrdiff_t h) const {
        if (empty()) {
            return {};
        }
        std::vector<size_t> res;
        bit_vector mask1(at(0).size(), true);
        mask1 += bit_vector(h, true);
        for (size_t i = 0; i < size(); ++i) {
            bit_vector const& v = at(i);
            if (v.all(mask1)) {
                res.push_back(i);
            }
        }
        return res;
    }

    std::vector<size_t> bit_matrix::get_g_p(ptrdiff_t h) const {
        if (empty()) {
            return {};
        }
        std::vector<size_t> res;
        bit_vector mask1(at(0).size(), false);
        mask1 += bit_vector(h, true);
        for (size_t i = 0; i < size(); ++i) {
            bit_vector const& v = at(i);
            if (v.all(mask1)) {
                res.push_back(i);
            }
        }
        return res;
    }

    std::vector<size_t> bit_matrix::get_g_f_s(ptrdiff_t f, ptrdiff_t s) const {
        if (empty() || s <= f) {
            return {};
        }
        std::vector<size_t> res;
        bit_vector mask1(at(0).size(), false);
        mask1 += bit_vector(s, true);
        bit_vector mask2(at(0).size(), true);
        mask2 += bit_vector(s, true);
        bit_vector mask3(at(0).size(), true);
        mask3 += bit_vector(f, true);

        for (size_t i = 0; i < size(); ++i) {
            bit_vector const& v = at(i);
            if (v.any(mask2) && v.any(mask1) && v.all(mask3)) {
                res.push_back(i);
            }
        }
        return res;
    }

    std::vector<size_t> bit_matrix::get_g_s_p(ptrdiff_t s, ptrdiff_t p) const {
        if (empty() || p <= s) {
            return {};
        }
        std::vector<size_t> res;
        bit_vector mask1(at(0).size(), false);
        mask1 += bit_vector(s, true);
        bit_vector mask2(at(0).size(), true);
        mask2 += bit_vector(s, true);
        bit_vector mask3(at(0).size(), false);
        mask3 += bit_vector(p, true);

        for (size_t i = 0; i < size(); ++i) {
            bit_vector const& v = at(i);
            if (v.any(mask2) && v.any(mask1) && v.all(mask3)) {
                res.push_back(i);
            }
        }
        return res;
    }

    bit_matrix bit_matrix::puncture(size_t x, size_t y) const {
        bit_matrix res;
        for (auto const& vect : *this) {
            res.emplace_back(vect.puncture(x, y));
        }
        return res;
    }



    bit_matrix bit_matrix::resolve_basis_gaussian() {
        bit_matrix res;
        for (size_t i = 0; i < size(); ++i) {
            size_t l = (*this)[i].leading();
            if (l != (*this)[i].size()) {
                res.push_back((*this)[i]);
                for (size_t j = i + 1; j < size(); ++j) {
                    if ((*this)[j][l]) {
                        (*this)[j] += (*this)[i];
                    }
                }
            }
        }
        return res;
    }

    bit_matrix bit_matrix::retrieve(std::vector<size_t> const& ind) const {
        bit_matrix res;
        for (size_t index : ind) {
            res.push_back(at(index));
        }
        return res;
    }

    bit_vector bit_matrix::get_and_multiply(bit_vector const& get) const {
        if (get.size() != size()) {
            throw std::logic_error("get and multiply: sizes are not matched");
        }
        bit_vector res(at(0).size(), true);
        for (size_t i = 0; i < size(); ++i) {
            if (get[i]) {
                res.multiply(at(i));
            }
        }
        return res;
    }

    size_t bit_matrix::get_c_tr_ctors_number(size_t h) const {
        size_t cnt = 0;
        bit_vector mask1(at(0).size(), false);
        mask1 += bit_vector(h, true);
        for (size_t i = 0; i < size(); ++i) {
            bit_vector const& v = at(i);
            if (v.all(mask1)) {
                ++cnt;
            }
        }
        return cnt;
    }

    void bit_matrix::make_tof() {
        for (size_t i = 0; i < size(); ++i) {
            size_t best_leading = 10e9;
            size_t best_number;
            for (size_t j = i; j < size(); ++j) {
                size_t leading = (*this)[j].leading();
                if (best_leading > leading) {
                    best_leading = leading;
                    best_number = j;
                }
            }
            if (best_leading == (*this)[i].size()) {
                throw std::logic_error("not a basis");
                break;
            }
            using std::swap;
            swap((*this)[i], (*this)[best_number]);
            for (size_t j = i + 1; j < size(); ++j) {
                if ((*this)[j][best_leading]) {
                    (*this)[j] += (*this)[i];
                }
            }
        }

        for (ptrdiff_t i = size() - 1; i > 0; --i) {
            for (ptrdiff_t j = i - 1; j >= 0; --j) {
                if ((*this)[j][(*this)[i].trailing()]) {
                    (*this)[j] += (*this)[i];
                }
            }
        }

    }

    std::string bit_matrix::to_string() const {
        std::string res;
        for (auto const& vect : *this) {
            res += vect.to_string();
            res += '\n';
        }
        return res;
    }


    lin_vector::lin_vector(size_t src, size_t dim) : std::vector<bool>(dim) {
        for (size_t i = 0; i < dim; ++i) {
            (*this)[i] = (src >> (dim - i - 1)) & 1;
        }
    }

    lin_vector& lin_vector::operator+=(lin_vector const& o) {
        if (size() != o.size()) {
            throw std::domain_error("cannot add vectors of different size");
        }
        for (size_t i = 0; i < size(); ++i) {
            (*this)[i] = (*this)[i] != o[i];
        }
        return *this;
    }

    lin_vector lin_vector::operator+(lin_vector const& o) const {
        lin_vector tmp(*this);
        tmp += o;
        return tmp;
    }

    int64_t lin_vector::operator*(lin_vector const& o) const {
        if (size() != o.size()) {
            throw std::domain_error("cannot multiply vectors of different size");
        }
        int64_t res = 0;
        for (size_t i = 0; i < size(); ++i) {
            res += ((*this)[i] && o[i]) ? 1 : 0;
        }
        return res;
    }

    lin_vector lin_vector::operator*(matrix const& o) const { // expects vector to be "horisontal"
        lin_vector res(o[0].size(), false);
        for (size_t i = 0; i < size(); ++i) {
            if (at(i)) {
                res += o[i];
            }
        }
        return res;
    }

    lin_vector& lin_vector::operator-() {
        for (size_t i = 0; i < size(); ++i) {
            at(i) = at(i) != true;
        }
        return *this;
    }

    lin_vector& lin_vector::multiply(lin_vector const& o) {
        if (size() != o.size()) {
            throw std::logic_error("cannot multiply vectors of different size");
        }
        for (size_t i = 0; i < size(); ++i) {
            at(i) = at(i) & o.at(i);
        }
        return *this;
    }

    void lin_vector::permutate(std::vector<size_t> const& permutation) {
        lin_vector res(size());

        for (size_t i = 0; i < size(); ++i) {
            res[i] = (*this)[permutation[i]];
        }
        *this = std::move(res);
    }

    matrix matrix::transpose() const {
        if (size() == 0) {
            return matrix();
        }
        matrix res((*this)[0].size(), lin_vector(size()));
        for (size_t i = 0; i < (*this)[0].size(); ++i) {
            for (size_t j = 0; j < size(); ++j) {
                res[i][j] = (*this)[j][i];
            }
        }
        return res;
    }

    lin_vector matrix::operator*(lin_vector const& o) const { // expects vector to be "vertical"
        if (size() == 0) {
            throw std::domain_error("cannot multiply with empty matrix");
        }
        if ((*this)[0].size() != o.size()) {
            throw std::domain_error("cannot multiply matrix and vector of different size");
        }
        lin_vector res(size(), false);
        for (size_t i = 0; i < size(); ++i) {
            for (size_t j = 0; j < (*this)[i].size(); ++j) {
                res[i] = res[i] != ((*this)[i][j] && o[j]);
            }
        }
        return res;
    }

    std::string lin_vector::to_string() const {
        std::string res;
        for (auto const& it : *this) {
            if (it) {
                res.push_back('1');
            }
            else {
                res.push_back('0');
            }
        }
        return res;
    }

    lin_vector& lin_vector::operator++() noexcept {
        for (ptrdiff_t i = size() - 1; i >= 0; --i) {
            if (!(*this)[i]) {
                (*this)[i] = true;
                return *this;
            }
            else {
                (*this)[i] = false;
            }
        }
        return *this;
    }

    uint64_t lin_vector::to_bit_mask() const {
        if (size() > 64) {
            throw std::domain_error("to big vector to convert to int");
        }
        uint64_t res = 0;
        for (size_t i = 0; i < size(); ++i) {
            res <<= 1;
            if ((*this)[i]) {
                ++res;
            }
        }
        return res;
    }


    size_t lin_vector::leading() const noexcept {
        size_t res = 0;
        for (; res < size(); ++res) {
            if (at(res)) {
                return res;
            }
        }
        return res;
    }


    size_t lin_vector::trailing() const noexcept {
        size_t res = size() - 1;
        for (; res >= 1; --res) {
            if (at(res)) {
                return res;
            }
        }
        return res;
    }

    bool lin_vector::all_zeros(size_t from, size_t to) const {
        bool flag = false;
        for (size_t i = from; i < to; ++i) {
            flag |= at(i);
        }
        return !flag;
    }

    lin_vector lin_vector::puncture(size_t x, size_t y) const {
        if (x > size() || y < x) {
            throw std::logic_error("bad bounds for puncturing");
        }
        lin_vector res;
        res.insert(res.begin(), begin() + x, begin() + y);
        return res;
    }

    lin_vector lin_vector::concat(lin_vector const& o) const {
        lin_vector res(*this);
        res.insert(res.end(), o.begin(), o.end());
        return res;
    }

    matrix matrix::operator*(matrix const& o) const {
        if (size() == 0 || o.size() == 0) {
            throw std::domain_error("cannot multiply with empty matrix");
        }
        if ((*this)[0].size() != o.size()) {
            throw std::domain_error("cannot multiply not matched matrices");
        }

        matrix res;
        res.resize(size());


        for (size_t i = 0; i < size(); ++i) {
            res[i].resize(o[0].size());
            for (size_t j = 0; j < o[0].size(); ++j) {
                res[i][j] = false;
                for (size_t k = 0; k < o.size(); ++k) {
                    res[i][j] = res[i][j] != ((*this)[i][k] && o[k][j]);
                }
            }
        }

        return res;
    }

    matrix matrix::inverse() const {

        matrix tmp = (*this);

        // prepare for gaussian method
        for (size_t i = 0; i < size(); ++i) {
            for (size_t j = 0; j < size(); ++j) {
                if (j == i) {
                    tmp[i].push_back(true);
                }
                else {
                    tmp[i].push_back(0);
                }
            }
        }

        // the basis on first size() vectors would be found, 
        // this matrix is k * 2k, at right was id matrix, 
        // so at the right there is inversed matrix found by gaussian method
        tmp.resolve_basis();
        matrix res(size());

        for (size_t i = 0; i < size(); ++i) {
            res[i].insert(res[i].end(), tmp[i].begin() + size(), tmp[i].end());
        }

        return res;
    }

    void matrix::permutate(std::vector<size_t> const& permutation) {
        for (auto& vect : *this) {
            vect.permutate(permutation);
        }
    }

    matrix matrix::resolve_basis_gaussian() {
        matrix res;
        for (size_t i = 0; i < size(); ++i) {
            if ((*this)[i].leading() != (*this)[i].size()) {
                res.push_back((*this)[i]);
                for (size_t j = i + 1; j < size(); ++j) {
                    if ((*this)[j][(*this)[i].leading()]) {
                        (*this)[j] += (*this)[i];
                    }
                }
            }
        }
        return res;
    }


    // returns vector of basis pos and transformation matrix to old basis
    std::pair<std::vector<size_t>, matrix> matrix::resolve_basis() {
        matrix old(*this);
        std::set<size_t> used_rows;
        size_t cnt = 0;
        std::vector<size_t> basis_pos;
        matrix transformation(size(), lin_vector(size(), false));


        basis_pos.reserve(size());
        for (size_t i = 0; i < (*this)[0].size(); ++i) {
            bool is_null_vect = true;
            size_t not_null_pos = 0;
            for (size_t j = 0; j < size(); ++j) { // find first non null vector
                if ((*this)[j][i] && (used_rows.find(j) == used_rows.end())) {
                    is_null_vect = false;
                    not_null_pos = j;
                    break;
                }
            }
            if (is_null_vect) {
                continue;
            }
            for (size_t j = 0; j < size(); ++j) {
                transformation[j][cnt] = old[j][i];
            }

            ++cnt;

            basis_pos.push_back(i); // add basis vector
            used_rows.insert(not_null_pos);
            for (size_t j = 0; j < size(); ++j) { // subtract basis vector if coord corresponding to it is equal to 1

                if (j == not_null_pos) {
                    continue;
                }
                if ((*this)[j][i]) {
                    (*this)[j] += (*this)[not_null_pos];
                }
            }


            if (cnt >= size()) { // if found basis then break
                break;
            }

        }

        // make identity matrix on basis vectors by permutations
        for (size_t i = 0; i < basis_pos.size(); ++i) {
            for (size_t j = 0; j < size(); ++j) {
                if ((*this)[j][basis_pos[i]]) {
                    using std::swap;
                    swap((*this)[j], (*this)[i]);
                }
            }
        }


        return std::make_pair(std::move(basis_pos), std::move(transformation));
    }

    void matrix::make_tof() {
        for (size_t i = 0; i < size(); ++i) {
            size_t best_leading = 10e9;
            size_t best_number;
            for (size_t j = i; j < size(); ++j) {
                size_t leading = (*this)[j].leading();
                if (best_leading > leading) {
                    best_leading = leading;
                    best_number = j;
                }
            }
            if (best_leading == (*this)[i].size()) {
                throw std::logic_error("not a basis");
                break;
            }
            using std::swap;
            swap((*this)[i], (*this)[best_number]);
            for (size_t j = i + 1; j < size(); ++j) {
                if ((*this)[j][best_leading]) {
                    (*this)[j] += (*this)[i];
                }
            }
        }

        for (ptrdiff_t i = size() - 1; i > 0; --i) {
            for (ptrdiff_t j = i - 1; j >= 0; --j) {
                if ((*this)[j][(*this)[i].trailing()]) {
                    (*this)[j] += (*this)[i];
                }
            }
        }

        if (!is_tof()) {
            throw std::logic_error("cant be not in tof");
        }

    }

    matrix matrix::retrieve(std::vector<size_t> const& ind) const {
        matrix res;
        for (size_t index : ind) {
            res.push_back(at(index));
        }
        return res;
    }

    std::vector<size_t> matrix::get_g_s(ptrdiff_t h) const {
        std::vector<size_t> res;
        for (size_t i = 0; i < size(); ++i) {
            if (h > at(i).leading() && h <= at(i).trailing()) {
                res.push_back(i);
            }
        }
        return res;
    }

    std::vector<size_t> matrix::get_g_f(ptrdiff_t h) const {
        std::vector<size_t> res;
        for (size_t i = 0; i < size(); ++i) {
            if (h <= at(i).leading()) {
                res.push_back(i);
            }
        }
        return res;
    }

    std::vector<size_t> matrix::get_g_p(ptrdiff_t h) const {
        std::vector<size_t> res;
        for (size_t i = 0; i < size(); ++i) {
            if (h > at(i).trailing()) {
                res.push_back(i);
            }
        }
        return res;
    }

    std::vector<size_t> matrix::get_g_f_s(ptrdiff_t f, ptrdiff_t s) const {
        if (s <= f) {
            return {};
        }
        auto indeces = get_g_f(f);
        auto snd = retrieve(get_g_f(f)).get_g_s(s);
        std::vector<size_t> res;
        for (size_t i = 0; i < snd.size(); ++i) {
            res.push_back(indeces[snd[i]]);
        }
        return res;
    }

    std::vector<size_t> matrix::get_g_s_p(ptrdiff_t f, ptrdiff_t s) const {
        if (s <= f) {
            return {};
        }
        auto indeces = get_g_s(f);
        auto snd = retrieve(get_g_s(f)).get_g_p(s);
        std::vector<size_t> res;
        for (size_t i = 0; i < snd.size(); ++i) {
            res.push_back(indeces[snd[i]]);
        }
        return res;
    }

    matrix matrix::puncture(size_t x, size_t y) const {
        matrix res;
        for (auto const& vect : *this) {
            res.emplace_back(vect.puncture(x, y));
        }
        return res;
    }

    std::string matrix::to_string() const {
        std::string res;
        for (auto const& vect : *this) {
            res += vect.to_string();
            res += '\n';
        }
        return res;
    }

    matrix matrix::operator+(matrix const& o) const {
        matrix mt(*this);
        for (size_t i = 0; i < size(); ++i) {
            mt[i] += o[i];
        }
        return mt;
    }

    lin_vector matrix::get_and_multiply(const lin_vector& get) const {
        if (get.size() != size()) {
            throw std::logic_error("get and multiply: sizes are not matched");
        }
        lin_vector res(at(0).size(), true);
        for (size_t i = 0; i < size(); ++i) {
            if (get.at(i)) {
                res.multiply(at(i));
            }
        }
        return res;
    }

    size_t matrix::get_c_tr_ctors_number(size_t t) const {
        size_t tr_ctors = 0;
        for (auto const& v : *this) {
            if (v.trailing() < t) {
                ++tr_ctors;
            }
        }
        return tr_ctors;
    }

    bool matrix::is_tof() const noexcept {
        ptrdiff_t prev = -1;

        for (auto const& v : *this) {
            if (((ptrdiff_t)v.leading()) <= prev) {
                return false;
            }
            else {
                prev = v.leading();
            }
        }

        std::set<size_t> last;

        for (auto const& v : *this) {
            if (last.find(v.trailing()) != last.end()) {
                return false;
            }
            else {
                last.insert(v.trailing());
            }
        }

        return true;

    }
}
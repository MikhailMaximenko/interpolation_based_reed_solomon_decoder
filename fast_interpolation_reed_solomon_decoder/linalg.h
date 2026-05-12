#pragma once

#include "socow_vector.h"

#include <cstddef>
#include <cstdint>
#include <functional>
#include <set>
#include <string>
#include <vector>

namespace linalg {

    struct matrix;
    struct bit_matrix;


    struct bit_vector {
        static constexpr int UNDERLIING_TYPE_SIZE = 8 * sizeof(uint64_t);
        socow_vector<uint64_t, 1> _storage;
        size_t _sz;

        bit_vector();
        bit_vector(size_t, bool v = false);
        bit_vector(size_t, size_t);
        explicit bit_vector(std::vector<bool> const&);


        void push_back(bool);
        size_t size() const noexcept;
        void resize(size_t);
        bool empty() const noexcept;
        bool operator[](size_t) const;
        void set(size_t, bool);

        void cp(bit_vector const&);


        bit_vector& operator+=(bit_vector const&);
        bit_vector operator+(bit_vector const&) const;

        bit_vector& add_with_offset(bit_vector const&, size_t);

        bit_vector operator-() const;

        bit_vector& multiply(bit_vector const&);

        void permutate(std::vector<size_t> const&);
        bit_vector operator*(bit_matrix const&) const;


        bit_vector& operator++() noexcept;

        bool any(bit_vector const& mask) const;
        bool all(bit_vector const& mask) const;

        size_t leading() const noexcept;
        size_t trailing() const noexcept;

        bool all_zeros(size_t, size_t) const;

        bit_vector puncture(size_t, size_t) const;
        bit_vector concat(bit_vector const&) const;
        bit_vector& concat(bit_vector const&);

        std::string to_string() const;

        bool operator==(bit_vector const& o) const noexcept;
        bool operator!=(bit_vector const& o) const noexcept;
        bool operator<(bit_vector const& o) const noexcept;

        uint64_t to_bit_mask() const noexcept;

    };

    struct bit_matrix : std::vector<bit_vector> {
        using std::vector<bit_vector>::vector;

        std::vector<size_t> get_g_s(ptrdiff_t) const;

        std::vector<size_t> get_g_f(ptrdiff_t) const;

        std::vector<size_t> get_g_p(ptrdiff_t) const;

        std::vector<size_t> get_g_f_s(ptrdiff_t, ptrdiff_t) const;

        std::vector<size_t> get_g_s_p(ptrdiff_t, ptrdiff_t) const;

        bit_matrix puncture(size_t, size_t) const;

        bit_matrix resolve_basis_gaussian();

        bit_matrix retrieve(std::vector<size_t> const&) const;

        bit_vector get_and_multiply(bit_vector const&) const;

        size_t get_c_tr_ctors_number(size_t) const;

        bool is_tof() const noexcept;

        void make_tof();

        std::string to_string() const;
    };

    struct lin_vector : std::vector<bool> {
        using std::vector<bool>::vector;

        lin_vector() = default;
        lin_vector(lin_vector const&) = default;
        lin_vector(lin_vector&&) = default;

        lin_vector& operator=(lin_vector const&) = default;
        lin_vector& operator=(lin_vector&&) = default;

        ~lin_vector() = default;


        lin_vector(size_t, size_t);

        lin_vector& operator+=(lin_vector const&);
        lin_vector operator+(lin_vector const&) const;
        int64_t operator*(lin_vector const&) const;

        lin_vector& operator-();

        lin_vector& multiply(lin_vector const&);

        void permutate(std::vector<size_t> const&);
        lin_vector operator*(matrix const&) const;


        lin_vector& operator++() noexcept;

        size_t leading() const noexcept;
        size_t trailing() const noexcept;

        bool all_zeros(size_t, size_t) const;

        lin_vector puncture(size_t, size_t) const;
        lin_vector concat(lin_vector const&) const;

        uint64_t to_bit_mask() const;
        std::string to_string() const;

    };

    struct matrix : std::vector<lin_vector> {
        using std::vector<lin_vector>::vector;

        matrix transpose() const;

        lin_vector operator*(lin_vector const&) const;
        matrix operator*(matrix const&) const;

        matrix inverse() const;

        void permutate(std::vector<size_t> const&);

        void make_tof();

        std::vector<size_t> get_g_s(ptrdiff_t) const;

        std::vector<size_t> get_g_f(ptrdiff_t) const;

        std::vector<size_t> get_g_p(ptrdiff_t) const;

        std::vector<size_t> get_g_f_s(ptrdiff_t, ptrdiff_t) const;

        std::vector<size_t> get_g_s_p(ptrdiff_t, ptrdiff_t) const;

        matrix puncture(size_t, size_t) const;

        matrix resolve_basis_gaussian();

        matrix retrieve(std::vector<size_t> const&) const;

        std::pair<std::vector<size_t>, matrix> resolve_basis();

        std::string to_string() const;

        matrix operator+(matrix const&) const;

        lin_vector get_and_multiply(lin_vector const&) const;

        size_t get_c_tr_ctors_number(size_t) const;

        bool is_tof() const noexcept;
    };

} // namespace linalg

template<>
struct std::hash<linalg::lin_vector> {
    std::hash<std::vector<bool>> hasher = std::hash<std::vector<bool>>();
    size_t operator()(const linalg::lin_vector& lv) const {
        return hasher(static_cast<std::vector<bool>const&>(lv));
    }
};

template <class T>
inline void hash_combine(std::size_t& seed, const T& v) {
    seed ^= std::hash<T>{}(v)+0x9e3779b9 + (seed << 6) + (seed >> 2);
}

template <>
struct std::hash<linalg::bit_vector> {
    std::size_t operator()(const linalg::bit_vector& v) const {
        std::size_t seed = v.size();
        for (uint64_t i : v._storage) {
            hash_combine(seed, i);
        }
        return seed;
    }
};
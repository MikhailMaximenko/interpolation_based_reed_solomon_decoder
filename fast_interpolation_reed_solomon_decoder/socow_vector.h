#pragma once
#include <algorithm>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <memory>
#include <utility>

template <typename T, size_t SMALL_SIZE>
class socow_vector {
public:
    using value_type = T;

    using reference = T&;
    using const_reference = const T&;

    using pointer = T*;
    using const_pointer = const T*;

    using iterator = pointer;
    using const_iterator = const_pointer;

private:
    struct dynamic_buffer {
        friend socow_vector;

        dynamic_buffer() : _cap(0), _refs(1) {}

    private:
        size_t _cap;
        int32_t _refs;
        T _data[0];
    };

    size_t _size;
    bool _small;

    union {
        T _static_data[SMALL_SIZE];
        dynamic_buffer* _dynamic_data;
    };

public:
    socow_vector() : _size(0), _small(true) {}

    socow_vector(const socow_vector& other) : _size(other._size), _small(other._small) {
        if (_small) {
            std::uninitialized_copy(other.begin(), other.end(), _static_data);
        }
        else {
            share(other._dynamic_data);
        }
    }

    socow_vector(size_t sz, const T& v) : _size(sz), _small(sz <= SMALL_SIZE) {
        if (!_small) {
            _dynamic_data = new (operator new(sizeof(dynamic_buffer) + sz * sizeof(T))) dynamic_buffer;
            for (size_t i = 0; i < sz; ++i) {
                _dynamic_data->_data[i] = v;
            }
        }
        else {
            for (size_t i = 0; i < sz; ++i) {
                _static_data[i] = v;
            }
        }
    }

    void swap(socow_vector& other) {
        if (this == &other) {
            return;
        }
        if (!_small && !other._small) {
            std::swap(_dynamic_data, other._dynamic_data);
        }
        else if (_small && !other._small) {
            other.swap(*this);
            return;
        }
        else if (!_small && other._small) {
            dynamic_buffer* tmp = _dynamic_data;
            try {
                std::uninitialized_copy_n(std::as_const(other).begin(), other.size(), _static_data);
            }
            catch (...) {
                _dynamic_data = tmp;
                throw;
            }
            std::destroy_n(other._static_data, other._size);
            other._dynamic_data = tmp;
        }
        else {
            if (_size < other._size) {
                other.swap(*this);
                return;
            }
            std::uninitialized_copy(_static_data + other._size, _static_data + _size, other._static_data + other._size);
            try {
                std::swap_ranges(_static_data, _static_data + other._size, other._static_data);
            }
            catch (...) {
                std::destroy(other._static_data + other._size, other._static_data + _size);
                throw;
            }
            std::destroy(_static_data + other._size, _static_data + _size);
        }
        std::swap(_size, other._size);
        std::swap(_small, other._small);
    }

    socow_vector& operator=(const socow_vector& other) {
        if (this == &other) {
            return *this;
        }
        if (other._small) {
            if (_small) {
                size_t min_size = std::min(_size, other._size);
                socow_vector temp(std::as_const(*this).begin(), std::as_const(*this).begin() + min_size, min_size);
                std::destroy(_static_data, _static_data + min_size);
                try {
                    std::uninitialized_copy(other._static_data, other._static_data + other._size, _static_data);
                    if (_size > other._size) {
                        std::destroy(_static_data + other._size, _static_data + _size);
                    }
                }
                catch (...) {
                    std::uninitialized_copy(temp.data(), temp.end(), _static_data);
                    throw;
                }
            }
            else {
                dynamic_buffer* tmp = _dynamic_data;
                try {
                    std::uninitialized_copy_n(other._static_data, other._size, _static_data);
                }
                catch (...) {
                    _dynamic_data = tmp;
                    throw;
                }
                unshare_dynamic_data(tmp);
            }
        }
        else {
            unshare();
            share(other._dynamic_data);
        }
        _size = other._size;
        _small = other._small;
        return *this;
    }

    void pop_back() {
        if (!_small) {
            if (_dynamic_data->_refs != 1) {
                socow_vector tmp(_dynamic_data->_data, _dynamic_data->_data + _size - 1, capacity());
                unshare();
                --_size;
                share(tmp._dynamic_data);
            }
            else {
                _dynamic_data->_data[--_size].~T();
            }

        }
        else {
            _static_data[--_size].~T();
        }
    }

    void push_back(const T& value) {
        insert(std::as_const(*this).end(), value);
    }

    bool empty() const noexcept {
        return _size == 0;
    }

    size_t capacity() const noexcept {
        return (_small) ? SMALL_SIZE : _dynamic_data->_cap;
    }

    void reserve(size_t new_cap) {
        if (new_cap < _size) {
            return;
        }
        if (!_small && _dynamic_data->_refs == 1 && capacity() >= new_cap) {
            return;
        }
        if (new_cap <= SMALL_SIZE) {
            shrink_to_fit();
            return;
        }
        socow_vector tmp(std::as_const(*this).begin(), std::as_const(*this).end(), new_cap);
        unshare();
        share(tmp._dynamic_data);
        _small = false;
    }

    void resize(size_t new_size, T const& v = T()) {
        if (_size >= new_size) {
            return;
        }
        reserve(new_size);
        for (size_t i = _size; i < new_size; ++i) {
            (*this)[i] = v;
        }
        _size = new_size;

    }

    void shrink_to_fit() {
        if (_small || _size == capacity()) {
            return;
        }
        if (_size <= SMALL_SIZE) {
            auto tmp = _dynamic_data;
            try {
                std::uninitialized_copy(std::as_const(*this).begin(), std::as_const(*this).end(), _static_data);
                unshare_dynamic_data(tmp);
                _small = true;
            }
            catch (...) {
                _dynamic_data = tmp;
                throw;
            }
        }
        else {
            socow_vector tmp(std::as_const(*this).begin(), std::as_const(*this).end(), _size);
            std::swap(_dynamic_data, tmp._dynamic_data);
        }
    }

    pointer data() {
        if (_small) {
            return _static_data;
        }
        copy_on_write();
        return _dynamic_data->_data;
    }

    const_pointer data() const noexcept {
        return _small ? _static_data : _dynamic_data->_data;
    }

    reference front() {
        return *data();
    }

    const_reference front() const noexcept {
        return *std::as_const(*this).data();
    }

    reference back() {
        return *(data() + _size - 1);
    }

    const_reference back() const noexcept {
        return *(std::as_const(*this).data() + _size - 1);
    }

    iterator begin() {
        return data();
    }

    const_iterator begin() const noexcept {
        return std::as_const(*this).data();
    }

    iterator end() {
        return data() + _size;
    }

    const_iterator end() const noexcept {
        return std::as_const(*this).data() + _size;
    }

    reference operator[](size_t ind) {
        return *(data() + ind);
    }

    const_reference operator[](size_t ind) const noexcept {
        return *(std::as_const(*this).data() + ind);
    }

    iterator insert(const_iterator pos, const T& value) {
        ptrdiff_t offset = pos - std::as_const(*this).data();
        if (_size == capacity() || (!_small && _dynamic_data->_refs != 1)) {
            int32_t cnst = (_size == capacity() ? 2 : 1);
            socow_vector tmp(std::as_const(*this).begin(), pos, cnst * capacity());
            tmp.push_back(value);
            std::uninitialized_copy(std::as_const(*this).begin() + offset, std::as_const(*this).end(),
                tmp.data() + offset + 1);
            unshare();
            share(tmp._dynamic_data);
            _small = false;
        }
        else {
            new (end()) T(value);
            for (int64_t i = _size; i > offset; --i) {
                try {
                    std::swap(data()[i], data()[i - 1]);
                }
                catch (...) {
                    for (int64_t j = i; j < _size; ++j) {
                        std::swap(data()[j], data()[j + 1]);
                    }
                    data()[_size].~T();
                    throw;
                }
            }
        }
        ++_size;
        return begin() + offset;
    }

    iterator erase(const_iterator pos) {
        return erase(pos, pos + 1);
    }

    iterator erase(const_iterator first, const_iterator last) {
        iterator data_begin = (_small ? _static_data : _dynamic_data->_data);
        ptrdiff_t offset1 = first - data_begin;
        ptrdiff_t offset2 = last - data_begin;
        if (!_small && _dynamic_data->_refs != 1) {
            socow_vector tmp(std::as_const(*this).begin(), first, capacity());
            std::uninitialized_copy(last, std::as_const(*this).end(), tmp.data() + offset1);
            unshare();
            share(tmp._dynamic_data);
            _size = _size - offset2 + offset1;
        }
        else {
            iterator it1 = begin() + offset1;
            iterator it2 = begin() + offset2;
            if (it1 != it2) {
                while (it2 != end()) {
                    std::iter_swap(it1, it2);
                    ++it1;
                    ++it2;
                }
                for (ptrdiff_t i = 0; i < last - first; ++i) {
                    pop_back();
                }
            }
        }
        return begin() + offset1;
    }

    size_t size() const noexcept {
        return _size;
    }

    void clear() {
        if (_small) {
            std::destroy(_static_data, _static_data + _size);
        }
        else {
            if (_dynamic_data->_refs == 1) {
                std::destroy_n(_dynamic_data->_data, _size);
            }
            else {
                unshare();
                _small = true;
            }
        }
        _size = 0;
    }

    ~socow_vector() {
        unshare();
    }

    bool operator==(socow_vector const& o) const& {
        if (_size != o._size) {
            return false;
        }
        auto it1 = begin();
        auto it2 = o.begin();
        while (it1 != end() && it2 != end()) {
            if (*it1 != *it2) {
                return false;
            }
            ++it1;
            ++it2;
        }
        return true;
    }

    bool operator!=(socow_vector const& o) const& {
        return !(*this == o);
    }

private:
    void copy_on_write() {
        if (_small || _dynamic_data->_refs == 1) {
            return;
        }
        auto tmp = copy_dynamic_buffer(_dynamic_data->_data, _size);
        --_dynamic_data->_refs;
        _dynamic_data = tmp;
    }

    void unshare() noexcept {
        if (!_small) {
            unshare_dynamic_data(_dynamic_data);
        }
        else {
            std::destroy_n(_static_data, _size);
        }
    }

    socow_vector(const_iterator begin, const_iterator end, size_t cap) : _size(end - begin), _small(cap <= SMALL_SIZE) {
        if (!_small) {
            _dynamic_data = new (operator new(sizeof(dynamic_buffer) + cap * sizeof(T))) dynamic_buffer;
            _dynamic_data->_cap = cap;
            _dynamic_data->_refs = 1;
            try {
                std::uninitialized_copy(begin, end, _dynamic_data->_data);
            }
            catch (...) {
                operator delete(_dynamic_data);
                throw;
            }
            return;
        }
        std::uninitialized_copy(begin, end, _static_data);
    }

    void unshare_dynamic_data(dynamic_buffer* src) {
        if (src->_refs == 1) {
            std::destroy_n(src->_data, _size);
            operator delete(src);
        }
        else {
            --src->_refs;
        }
    }

    socow_vector(const socow_vector& other, size_t cap) : _size(other._size), _small(cap <= SMALL_SIZE) {
        assert(_size <= cap);

        if (_small) {
            std::uninitialized_copy(other.begin(), other.end(), _static_data);
            return;
        }
        if (cap == other._dynamic_data->_cap) {
            share(other._dynamic_data);
            return;
        }
        _dynamic_data = copy_dynamic_buffer(std::as_const(other).data(), cap, _size);
    }

    void share(dynamic_buffer* src) {
        _dynamic_data = src;
        ++src->_refs;
    }

    dynamic_buffer* copy_dynamic_buffer(const T* buf, size_t cap, size_t size) {
        auto tmp = new (operator new(sizeof(dynamic_buffer) + cap * sizeof(T))) dynamic_buffer;
        try {
            std::uninitialized_copy(buf, buf + size, tmp->_data);
        }
        catch (...) {
            operator delete(tmp);
            throw;
        }
        tmp->_cap = cap;
        return tmp;
    }

    dynamic_buffer* copy_dynamic_buffer(const T* buf, size_t cap) {
        return copy_dynamic_buffer(buf, cap, _size);
    }
};
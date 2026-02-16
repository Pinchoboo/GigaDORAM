//! should only be included internally by utils

#pragma once
#include "globals.h"

typedef uint32_t uint;

//! all of the functions below should be checked for correctness
template <typename T> constexpr inline T get_all_ones_of_type() // ? the constexpr is highly questionable
{
    T all_ones{};
    auto *bytes = reinterpret_cast<uint8_t *>(&all_ones);
    for (size_t i = 0; i < sizeof(T); i++) {
        bytes[i] = 0xFF;
    }
    return all_ones;
}

template <typename T> constexpr inline T get_all_zeros_of_type()
{
    return T{};
}

template <typename Big, typename Small> inline Big _get_last_typesize_nuker() // was constexpr
{
    Big big_all_ones = get_all_ones_of_type<Big>();
    Small *tmp = (Small *)&big_all_ones;
    tmp[sizeof(Big) / sizeof(Small) - 1] = Small{};
    return *((Big *)tmp);
}

template <typename q_type, typename x_type> constexpr inline q_type get_lowest_x_type_many_bits_are_ones()
{
    q_type ones_arr{};
    auto *bytes = reinterpret_cast<uint8_t *>(&ones_arr);
    for (size_t i = 0; i < sizeof(x_type); i++) {
        bytes[i] = 0xFF;
    }
    return ones_arr;
}

//! this only works for uint while the rest of the class is more general, it would be nice to generalize return val
template <typename T> uint inline __cht_for_oht_hash(const T *el, bool col, uint num_bits)
{
    assert(sizeof(T) * 8 > num_bits && num_bits < 32 && "these are the ranges this function works for");
    uint chunk = ((uint *)el)[col];
    uint mask = ~(UINT_MAX << num_bits);
    uint range_bit = ((uint)col) << num_bits;
    return (chunk & mask) | range_bit;
}

template <typename T> uint inline get_last_uint_bits(const T &el)
{
    assert(sizeof(T) % sizeof(uint) == 0);
    uint *uint_arr = (uint *)&el;
    return uint_arr[sizeof(T) / sizeof(uint) - 1];
}

template <typename T> T get_one_hot_mask (uint set_bit_index) {
    assert(set_bit_index < (8 * sizeof(T)));
    T mask{};
    ((char*)(&mask))[set_bit_index/8] |= (1 << (set_bit_index % 8));
    return mask;
}

template <typename T> T get_all_zero_except_nth_from_highest(uint n) {
    assert(n < (8 * sizeof(T)));
    return get_one_hot_mask<T>(8 * sizeof(T) - 1 - n);
}

template <typename T> T get_all_ones_rightshifted_by(uint shift) {
    assert(shift <= 8 * sizeof(T));
    T mask = get_all_ones_of_type<T>();
    for (uint i = 0; i < shift; i++) {
        mask ^= get_all_zero_except_nth_from_highest<T>(i);
    }
    return mask;
}

template <size_t YBytes>
inline emp::y_type_t<YBytes> get_all_ones_of_y_type() {
    emp::y_type_t<YBytes> all_ones{};
    const uint n = emp::y_type_t<YBytes>::active_bytes();
    for (uint i = 0; i < n; i++) {
        all_ones.data[i] = 0xFF;
    }
    return all_ones;
}

template <size_t YBytes>
inline emp::y_type_t<YBytes> get_one_hot_mask_y_type(uint set_bit_index) {
    assert(set_bit_index < emp::y_type_t<YBytes>::active_bits());
    emp::y_type_t<YBytes> mask{};
    mask.data[set_bit_index / 8] |= static_cast<uint8_t>(1u << (set_bit_index % 8));
    return mask;
}

template <size_t YBytes>
inline emp::y_type_t<YBytes> get_all_zero_except_nth_from_highest_y_type(uint n) {
    assert(n < emp::y_type_t<YBytes>::active_bits());
    return get_one_hot_mask_y_type<YBytes>(emp::y_type_t<YBytes>::active_bits() - 1 - n);
}

template <size_t YBytes>
inline emp::y_type_t<YBytes> get_all_ones_rightshifted_by_y_type(uint shift) {
    assert(shift <= emp::y_type_t<YBytes>::active_bits());
    emp::y_type_t<YBytes> mask = get_all_ones_of_y_type<YBytes>();
    for (uint i = 0; i < shift; i++) {
        mask ^= get_all_zero_except_nth_from_highest_y_type<YBytes>(i);
    }
    return mask;
}

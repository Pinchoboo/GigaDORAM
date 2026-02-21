#pragma once

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <stdexcept>
#include <vector>

#include "runtime_width.h"

namespace emp {
namespace runtime_y_ops {

inline void set_alibi_bit(uint8_t* bytes, uint32_t level, const runtime_width::RuntimeWidthSpec& spec);
inline bool get_bit_from_y_window(const uint8_t* bytes, uint32_t bit_index, const runtime_width::RuntimeWidthSpec& spec);

inline size_t checked_y_window_offset_bytes(
    uint64_t element_index,
    const runtime_width::RuntimeWidthSpec& spec,
    const char* ctx) {
    return static_cast<size_t>(runtime_width::checked_u64_mul(
        element_index,
        static_cast<uint64_t>(spec.y_stride_bytes),
        ctx
    ));
}

inline void clear_y_window(uint8_t* dst, const runtime_width::RuntimeWidthSpec& spec) {
    std::memset(dst, 0, spec.y_stride_bytes);
}

inline void copy_payload_into_y_window(
    uint8_t* dst,
    const uint8_t* src_payload,
    const runtime_width::RuntimeWidthSpec& spec) {
    clear_y_window(dst, spec);
    std::memcpy(dst, src_payload, spec.payload_bytes);
}

inline void copy_payload_from_y_window(
    uint8_t* dst_payload,
    const uint8_t* src,
    const runtime_width::RuntimeWidthSpec& spec) {
    std::memcpy(dst_payload, src, spec.payload_bytes);
}

inline void xor_y_window(uint8_t* dst, const uint8_t* src, const runtime_width::RuntimeWidthSpec& spec) {
    for (uint32_t i = 0; i < spec.y_stride_bytes; ++i) {
        dst[i] ^= src[i];
    }
}

inline void keep_payload_only(uint8_t* bytes, const runtime_width::RuntimeWidthSpec& spec) {
    if (spec.payload_bytes < spec.y_stride_bytes) {
        std::memset(bytes + spec.payload_bytes, 0, spec.y_stride_bytes - spec.payload_bytes);
    }
}

inline void keep_payload_and_alibi(uint8_t* bytes, const runtime_width::RuntimeWidthSpec& spec) {
    if (spec.alibi_bits == 0) {
        keep_payload_only(bytes, spec);
        return;
    }
    std::vector<uint8_t> keep_alibi(spec.alibi_bits, 0);
    for (uint32_t level = 0; level < spec.alibi_bits; ++level) {
        const uint32_t bit_index = spec.total_bits - 1 - level;
        keep_alibi[level] = get_bit_from_y_window(bytes, bit_index, spec) ? 1U : 0U;
    }
    keep_payload_only(bytes, spec);
    for (uint32_t level = 0; level < spec.alibi_bits; ++level) {
        if (keep_alibi[level] != 0) {
            set_alibi_bit(bytes, level, spec);
        }
    }
}

inline void set_alibi_bit(uint8_t* bytes, uint32_t level, const runtime_width::RuntimeWidthSpec& spec) {
    if (level >= spec.alibi_bits) {
        throw std::out_of_range("set_alibi_bit level out of range");
    }
    const uint32_t bit_index = spec.total_bits - 1 - level;
    const uint32_t byte_index = bit_index / 8;
    const uint8_t bit_mask = static_cast<uint8_t>(1u << (bit_index % 8));
    if (byte_index >= spec.y_stride_bytes) {
        throw std::out_of_range("set_alibi_bit byte index out of range");
    }
    bytes[byte_index] |= bit_mask;
}

inline bool get_bit_from_y_window(const uint8_t* bytes, uint32_t bit_index, const runtime_width::RuntimeWidthSpec& spec) {
    const uint32_t byte_index = bit_index / 8;
    if (byte_index >= spec.y_stride_bytes) {
        throw std::out_of_range("get_bit_from_y_window bit index out of range");
    }
    return ((bytes[byte_index] >> (bit_index % 8)) & 1u) != 0;
}

} // namespace runtime_y_ops
} // namespace emp

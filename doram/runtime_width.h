#pragma once

#include <cstdint>
#include <stdexcept>

#include "globals.h"

namespace emp {
namespace runtime_width {

struct RuntimeWidthSpec {
    uint32_t payload_bits = 0;
    uint32_t payload_bytes = 0;
    uint32_t alibi_bits = 0;
    uint32_t total_bits = 0;
    uint32_t packed_elem_bytes = 0;
    uint32_t packed_elem_blocks = 0;
};

inline constexpr uint32_t packed_elem_bytes_for_payload_bytes(uint32_t payload_bytes) {
    return static_cast<uint32_t>(2 * sizeof(x_type) + payload_bytes);
}

inline constexpr uint32_t packed_elem_blocks_for_payload_bytes(uint32_t payload_bytes) {
    return static_cast<uint32_t>(
        (packed_elem_bytes_for_payload_bytes(payload_bytes) + sizeof(block) - 1) / sizeof(block)
    );
}

inline RuntimeWidthSpec make_runtime_width_spec(
    uint32_t payload_bits,
    uint32_t alibi_bits,
    uint32_t total_bits) {
    if (payload_bits == 0 || (payload_bits % 8) != 0) {
        throw std::invalid_argument("payload_bits must be > 0 and divisible by 8");
    }
    if (alibi_bits >= total_bits) {
        throw std::invalid_argument("alibi_bits must be less than total_bits");
    }
    if (payload_bits > total_bits || (payload_bits + alibi_bits) > total_bits) {
        throw std::invalid_argument("payload_bits + alibi_bits must be <= total_bits");
    }

    RuntimeWidthSpec spec;
    spec.payload_bits = payload_bits;
    spec.payload_bytes = payload_bits / 8;
    spec.alibi_bits = alibi_bits;
    spec.total_bits = total_bits;
    spec.packed_elem_bytes = packed_elem_bytes_for_payload_bytes(spec.payload_bytes);
    spec.packed_elem_blocks = packed_elem_blocks_for_payload_bytes(spec.payload_bytes);
    return spec;
}

static_assert(packed_elem_bytes_for_payload_bytes(8) == 16, "runtime packed bytes 8B mismatch");
static_assert(packed_elem_bytes_for_payload_bytes(16) == 24, "runtime packed bytes 16B mismatch");
static_assert(packed_elem_bytes_for_payload_bytes(32) == 40, "runtime packed bytes 32B mismatch");
static_assert(packed_elem_bytes_for_payload_bytes(64) == 72, "runtime packed bytes 64B mismatch");
static_assert(packed_elem_bytes_for_payload_bytes(128) == 136, "runtime packed bytes 128B mismatch");

static_assert(packed_elem_blocks_for_payload_bytes(8) == 1, "runtime packed blocks 8B mismatch");
static_assert(packed_elem_blocks_for_payload_bytes(16) == 2, "runtime packed blocks 16B mismatch");
static_assert(packed_elem_blocks_for_payload_bytes(32) == 3, "runtime packed blocks 32B mismatch");
static_assert(packed_elem_blocks_for_payload_bytes(64) == 5, "runtime packed blocks 64B mismatch");
static_assert(packed_elem_blocks_for_payload_bytes(128) == 9, "runtime packed blocks 128B mismatch");

} // namespace runtime_width
} // namespace emp

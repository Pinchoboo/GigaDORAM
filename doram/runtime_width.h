#pragma once

#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>

#include "globals.h"

namespace emp {
namespace runtime_width {

struct RuntimeWidthSpec {
    uint32_t payload_bits = 0;
    uint32_t payload_bytes = 0;
    uint32_t alibi_bits = 0;
    uint32_t total_bits = 0;
    uint32_t y_blocks = 0;
    uint32_t y_stride_bytes = 0;
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

inline constexpr uint32_t y_blocks_for_payload_bytes(uint32_t payload_bytes) {
    return static_cast<uint32_t>((payload_bytes + sizeof(block) - 1) / sizeof(block));
}

inline constexpr uint32_t y_stride_bytes_for_payload_bytes(uint32_t payload_bytes) {
    return static_cast<uint32_t>(y_blocks_for_payload_bytes(payload_bytes) * sizeof(block));
}

inline constexpr uint32_t y_blocks_for_total_bits(uint32_t total_bits) {
    return static_cast<uint32_t>(((total_bits + 7U) / 8U + sizeof(block) - 1) / sizeof(block));
}

inline constexpr uint32_t y_stride_bytes_for_total_bits(uint32_t total_bits) {
    return static_cast<uint32_t>(y_blocks_for_total_bits(total_bits) * sizeof(block));
}

inline uint32_t checked_u32_add(uint32_t a, uint32_t b, const char* ctx) {
    if (b > (std::numeric_limits<uint32_t>::max() - a)) {
        throw std::overflow_error(std::string("u32 add overflow in ") + ctx);
    }
    return a + b;
}

inline uint32_t checked_u32_mul(uint32_t a, uint32_t b, const char* ctx) {
    if (a == 0 || b == 0) {
        return 0;
    }
    if (a > (std::numeric_limits<uint32_t>::max() / b)) {
        throw std::overflow_error(std::string("u32 mul overflow in ") + ctx);
    }
    return a * b;
}

inline uint64_t checked_u64_mul(uint64_t a, uint64_t b, const char* ctx) {
    if (a == 0 || b == 0) {
        return 0;
    }
    if (a > (std::numeric_limits<uint64_t>::max() / b)) {
        throw std::overflow_error(std::string("u64 mul overflow in ") + ctx);
    }
    return a * b;
}

inline uint64_t checked_u64_add(uint64_t a, uint64_t b, const char* ctx) {
    if (b > (std::numeric_limits<uint64_t>::max() - a)) {
        throw std::overflow_error(std::string("u64 add overflow in ") + ctx);
    }
    return a + b;
}

inline uint64_t checked_total_share_bytes(
    uint64_t batch_size,
    uint32_t value_stride_bytes,
    uint32_t rep_share_count,
    const char* ctx) {
    const uint64_t one_batch = checked_u64_mul(
        static_cast<uint64_t>(value_stride_bytes),
        static_cast<uint64_t>(rep_share_count),
        ctx
    );
    return checked_u64_mul(batch_size, one_batch, ctx);
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
    spec.y_blocks = y_blocks_for_total_bits(spec.total_bits);
    spec.y_stride_bytes = y_stride_bytes_for_total_bits(spec.total_bits);
    spec.packed_elem_bytes = packed_elem_bytes_for_payload_bytes(spec.payload_bytes);
    spec.packed_elem_blocks = packed_elem_blocks_for_payload_bytes(spec.payload_bytes);
    if (spec.y_blocks == 0 || spec.y_stride_bytes == 0) {
        throw std::invalid_argument("runtime width produced zero y block/stride");
    }
    if (spec.payload_bytes > spec.y_stride_bytes) {
        throw std::invalid_argument("runtime width payload bytes exceed y stride bytes");
    }
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

static_assert(y_blocks_for_payload_bytes(1) == 1, "runtime y blocks 1B mismatch");
static_assert(y_blocks_for_payload_bytes(16) == 1, "runtime y blocks 16B mismatch");
static_assert(y_blocks_for_payload_bytes(17) == 2, "runtime y blocks 17B mismatch");

static_assert(y_stride_bytes_for_payload_bytes(1) == 16, "runtime y stride 1B mismatch");
static_assert(y_stride_bytes_for_payload_bytes(16) == 16, "runtime y stride 16B mismatch");
static_assert(y_stride_bytes_for_payload_bytes(17) == 32, "runtime y stride 17B mismatch");

} // namespace runtime_width
} // namespace emp

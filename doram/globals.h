#pragma once

#include "emp-tool/emp-tool.h"
#include "rep_net_io_channel.h"
#include <thread>
#include <map> 


using namespace std;

namespace emp
{

class AbandonIO : public IOChannel<AbandonIO>
{
  public:
    void flush()
    {
    }

    void send_data_internal(const void *data, int len)
    {
    }

    void recv_data_internal(void *data, int len)
    {
    }
};

class BristolFashion_array;

// party conventions:
// 1: Party 1
// 2: Party 2
// 3: Party 3
int party = 0; 

uint NUM_THREADS = 2;

// Each of the thread_unsafe:: resources is an array of length num_threads
// which has to be initialized in main. The 0th element is the default
// for code that will never run in a thread. Code that runs in threads
// should look up its thread's assigned resources.
// The ith thread_unsafe::prev resource is assumed to be connected to
// the previous party's ith thread_unsafe::next resource.

class SHRepArray;
class MalRepArray;

PRG **prev_prgs;
PRG **next_prgs;
PRG **private_prgs;
PRG **shared_prgs;

RepNetIO **prev_ios;
RepNetIO **next_ios;
SHRepArray **rep_execs;


namespace thread_unsafe {

PRG *prev_prg;
PRG *next_prg;
// does private_prg really need to be thread_unsafe?
PRG *private_prg;
PRG *shared_prg;

RepNetIO *prev_io;
RepNetIO *next_io;
SHRepArray *rep_exec;

}

// TIMING
fstream timing_file;
fstream special_debug_file;

// ignore setup
double time_total = 0;         
vector<double> time_total_builds;
double time_total_build_prf = 0;
double time_total_batcher = 0;
double time_total_deletes = 0; 
double time_total_queries = 0;
double time_total_query_prf = 0;
double time_total_query_stupid = 0;
double time_total_shuffles = 0;
double time_total_cht_build = 0;
double time_total_transpose = 0;

// The following timings are not thread safe and commented out for now
/*
double time_total_mand = 0;
double time_compute_start = 0;
double time_compute_finish = 0;
double time_main_loop = 0;
*/

// double time_total_network = 0; // defined in rep_net_io to avoid include problem, can still be accessed globaly

double time_doram_constructor = 0; //? what is this for?

typedef unsigned long long ull;

typedef uint32_t x_type;

// ---------------------------------------------------------------------------
// y_type: fixed compile-time width.
// Override Y_TYPE_BITS at compile time to set the maximum width for builds.
// ---------------------------------------------------------------------------
#ifndef Y_TYPE_BITS
#define Y_TYPE_BITS 64
#endif
#define Y_TYPE_MAX_BITS Y_TYPE_BITS
#define Y_TYPE_MAX_BYTES ((Y_TYPE_MAX_BITS + 7) / 8)

template <size_t FixedBytes>
struct y_type_t {
    static_assert(FixedBytes > 0, "y_type_t requires FixedBytes > 0");
    uint8_t data[FixedBytes];

    static inline uint32_t active_bytes() {
        return static_cast<uint32_t>(FixedBytes);
    }

    static inline uint32_t active_bits() {
        return static_cast<uint32_t>(8 * FixedBytes);
    }

    // Default: zero
    y_type_t() : data{} {}

    // Implicit conversion from uint64_t (little-endian)
    y_type_t(uint64_t val) : data{} {
        const size_t n = active_bytes();
        for (size_t i = 0; i < sizeof(val) && i < n; i++)
            data[i] = static_cast<uint8_t>(val >> (i * 8));
    }

    // Bitwise XOR
    y_type_t operator^(const y_type_t& o) const {
        y_type_t r;
        const size_t n = active_bytes();
        for (size_t i = 0; i < n; i++) r.data[i] = data[i] ^ o.data[i];
        return r;
    }
    y_type_t& operator^=(const y_type_t& o) {
        const size_t n = active_bytes();
        for (size_t i = 0; i < n; i++) data[i] ^= o.data[i];
        return *this;
    }

    // Bitwise OR
    y_type_t operator|(const y_type_t& o) const {
        y_type_t r;
        const size_t n = active_bytes();
        for (size_t i = 0; i < n; i++) r.data[i] = data[i] | o.data[i];
        return r;
    }
    y_type_t& operator|=(const y_type_t& o) {
        const size_t n = active_bytes();
        for (size_t i = 0; i < n; i++) data[i] |= o.data[i];
        return *this;
    }

    // Bitwise AND
    y_type_t operator&(const y_type_t& o) const {
        y_type_t r;
        const size_t n = active_bytes();
        for (size_t i = 0; i < n; i++) r.data[i] = data[i] & o.data[i];
        return r;
    }
    y_type_t& operator&=(const y_type_t& o) {
        const size_t n = active_bytes();
        for (size_t i = 0; i < n; i++) data[i] &= o.data[i];
        return *this;
    }

    // Comparison
    bool operator==(const y_type_t& o) const {
        const size_t n = active_bytes();
        return memcmp(data, o.data, n) == 0;
    }
    bool operator!=(const y_type_t& o) const { return !(*this == o); }

    // Extract the low 64 bits as uint64_t (for B2A/A2B which stay 64-bit)
    uint64_t to_u64() const {
        uint64_t v = 0;
        const size_t n = active_bytes();
        for (size_t i = 0; i < sizeof(v) && i < n; i++)
            v |= static_cast<uint64_t>(data[i]) << (i * 8);
        return v;
    }
};

using y_type = y_type_t<Y_TYPE_MAX_BYTES>;

static_assert(sizeof(y_type_t<1>) == 1, "y_type_t size mismatch");
static_assert(sizeof(y_type_t<8>) == 8, "y_type_t size mismatch");
static_assert(sizeof(y_type_t<16>) == 16, "y_type_t size mismatch");
static_assert(sizeof(y_type_t<32>) == 32, "y_type_t size mismatch");
static_assert(sizeof(y_type_t<64>) == 64, "y_type_t size mismatch");
static_assert(sizeof(y_type_t<128>) == 128, "y_type_t size mismatch");

template <size_t FixedBytes>
inline constexpr uint32_t y_type_bits_v = static_cast<uint32_t>(8 * FixedBytes);

template <size_t FixedBytes>
inline constexpr uint32_t y_type_bytes_v = static_cast<uint32_t>(FixedBytes);

template <typename YType>
inline constexpr uint32_t y_type_bits_of() {
    return static_cast<uint32_t>(8 * sizeof(YType));
}

template <typename YType>
inline constexpr uint32_t y_type_bytes_of() {
    return static_cast<uint32_t>(sizeof(YType));
}

// Packed element layout used by xy_if_xs_equal and compare_swap circuits:
//   [x_type (4B)] [x_type (4B)] [y_type (compile-time bytes)] [padding]
template <typename YType>
inline constexpr uint packed_xy_bytes_for() {
    return 2 * sizeof(x_type) + y_type_bytes_of<YType>();
}
template <typename YType>
inline constexpr uint blocks_per_packed_xy_for() {
    return (packed_xy_bytes_for<YType>() + sizeof(block) - 1) / sizeof(block);
}
template <typename YType>
inline constexpr uint packed_xy_stride_for() {
    return blocks_per_packed_xy_for<YType>() * sizeof(block);
}

static_assert(packed_xy_bytes_for<y_type_t<8>>() == 16, "packed_xy_bytes_for<8> mismatch");
static_assert(packed_xy_bytes_for<y_type_t<16>>() == 24, "packed_xy_bytes_for<16> mismatch");
static_assert(packed_xy_bytes_for<y_type_t<32>>() == 40, "packed_xy_bytes_for<32> mismatch");
static_assert(packed_xy_bytes_for<y_type_t<64>>() == 72, "packed_xy_bytes_for<64> mismatch");
static_assert(packed_xy_bytes_for<y_type_t<128>>() == 136, "packed_xy_bytes_for<128> mismatch");

static_assert(blocks_per_packed_xy_for<y_type_t<8>>() == 1, "blocks_per_packed_xy_for<8> mismatch");
static_assert(blocks_per_packed_xy_for<y_type_t<16>>() == 2, "blocks_per_packed_xy_for<16> mismatch");
static_assert(blocks_per_packed_xy_for<y_type_t<32>>() == 3, "blocks_per_packed_xy_for<32> mismatch");
static_assert(blocks_per_packed_xy_for<y_type_t<64>>() == 5, "blocks_per_packed_xy_for<64> mismatch");
static_assert(blocks_per_packed_xy_for<y_type_t<128>>() == 9, "blocks_per_packed_xy_for<128> mismatch");

static_assert(packed_xy_stride_for<y_type_t<8>>() == 16, "packed_xy_stride_for<8> mismatch");
static_assert(packed_xy_stride_for<y_type_t<16>>() == 32, "packed_xy_stride_for<16> mismatch");
static_assert(packed_xy_stride_for<y_type_t<32>>() == 48, "packed_xy_stride_for<32> mismatch");
static_assert(packed_xy_stride_for<y_type_t<64>>() == 80, "packed_xy_stride_for<64> mismatch");
static_assert(packed_xy_stride_for<y_type_t<128>>() == 144, "packed_xy_stride_for<128> mismatch");

// Stream output for debugging
template <size_t FixedBytes>
inline std::ostream& operator<<(std::ostream& os, const y_type_t<FixedBytes>& v) {
    os << "0x";
    for (int i = static_cast<int>(y_type_t<FixedBytes>::active_bytes()) - 1; i >= 0; i--) {
        char buf[3];
        snprintf(buf, sizeof(buf), "%02x", v.data[i]);
        os << buf;
    }
    return os;
}

BristolFashion_array *cht_lookup_circuit_file = nullptr;
BristolFashion_array *prf_circuit = nullptr;
BristolFashion_array *replace_if_dummy_circuit_file[32];
BristolFashion_array *dummy_check_circuit_file[32];
BristolFashion_array *b2a_circuit_file = nullptr;
BristolFashion_array *a2b_circuit_file = nullptr;

} // namespace emp

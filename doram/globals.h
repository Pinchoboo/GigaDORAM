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
// y_type: compile-time configurable wide data type for ORAM slots.
// Override Y_TYPE_BITS at compile time (e.g. -DY_TYPE_BITS=272) to widen.
// Default: 64 bits (backward-compatible with original uint64_t y_type).
// ---------------------------------------------------------------------------
#ifndef Y_TYPE_BITS
#define Y_TYPE_BITS 64
#endif
#define Y_TYPE_BYTES ((Y_TYPE_BITS + 7) / 8)

struct y_type {
    uint8_t data[Y_TYPE_BYTES];

    // Default: zero
    y_type() { memset(data, 0, Y_TYPE_BYTES); }

    // Implicit conversion from uint64_t (little-endian)
    y_type(uint64_t val) {
        memset(data, 0, Y_TYPE_BYTES);
        for (size_t i = 0; i < sizeof(val) && i < Y_TYPE_BYTES; i++)
            data[i] = static_cast<uint8_t>(val >> (i * 8));
    }

    // Bitwise XOR
    y_type operator^(const y_type& o) const {
        y_type r;
        for (int i = 0; i < Y_TYPE_BYTES; i++) r.data[i] = data[i] ^ o.data[i];
        return r;
    }
    y_type& operator^=(const y_type& o) {
        for (int i = 0; i < Y_TYPE_BYTES; i++) data[i] ^= o.data[i];
        return *this;
    }

    // Bitwise OR
    y_type operator|(const y_type& o) const {
        y_type r;
        for (int i = 0; i < Y_TYPE_BYTES; i++) r.data[i] = data[i] | o.data[i];
        return r;
    }
    y_type& operator|=(const y_type& o) {
        for (int i = 0; i < Y_TYPE_BYTES; i++) data[i] |= o.data[i];
        return *this;
    }

    // Bitwise AND
    y_type operator&(const y_type& o) const {
        y_type r;
        for (int i = 0; i < Y_TYPE_BYTES; i++) r.data[i] = data[i] & o.data[i];
        return r;
    }
    y_type& operator&=(const y_type& o) {
        for (int i = 0; i < Y_TYPE_BYTES; i++) data[i] &= o.data[i];
        return *this;
    }

    // Comparison
    bool operator==(const y_type& o) const { return memcmp(data, o.data, Y_TYPE_BYTES) == 0; }
    bool operator!=(const y_type& o) const { return !(*this == o); }

    // Extract the low 64 bits as uint64_t (for B2A/A2B which stay 64-bit)
    uint64_t to_u64() const {
        uint64_t v = 0;
        for (size_t i = 0; i < sizeof(v) && i < Y_TYPE_BYTES; i++)
            v |= static_cast<uint64_t>(data[i]) << (i * 8);
        return v;
    }
};

// Packed element layout used by xy_if_xs_equal and compare_swap circuits:
//   [x_type (4B)] [x_type (4B)] [y_type (Y_TYPE_BYTES)] [padding to block boundary]
// Number of bytes of payload:
constexpr uint PACKED_XY_BYTES = 2 * sizeof(x_type) + Y_TYPE_BYTES;
// Number of 128-bit blocks per packed element:
constexpr uint BLOCKS_PER_PACKED_XY = (PACKED_XY_BYTES + sizeof(block) - 1) / sizeof(block);
// Byte stride per packed element (aligned to block boundary):
constexpr uint PACKED_XY_STRIDE = BLOCKS_PER_PACKED_XY * sizeof(block);

// Stream output for debugging
inline std::ostream& operator<<(std::ostream& os, const y_type& v) {
    os << "0x";
    for (int i = Y_TYPE_BYTES - 1; i >= 0; i--) {
        char buf[3];
        snprintf(buf, sizeof(buf), "%02x", v.data[i]);
        os << buf;
    }
    return os;
}

BristolFashion_array *xy_if_xs_equal_circuit = nullptr;
BristolFashion_array *cht_lookup_circuit_file = nullptr;
BristolFashion_array *prf_circuit = nullptr;
BristolFashion_array *replace_if_dummy_circuit_file[32];
BristolFashion_array *dummy_check_circuit_file[32];
BristolFashion_array *compare_swap_circuit_file = nullptr;
BristolFashion_array *b2a_circuit_file = nullptr;
BristolFashion_array *a2b_circuit_file = nullptr;

} // namespace emp

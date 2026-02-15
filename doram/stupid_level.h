#pragma once

#include "globals.h"
#include "rep_array_unsliced.h"
#include "bristol_fashion_array.h"
#include "utils.h"

#include "emp-tool/emp-tool.h"

namespace emp
{

class StupidLevel
{
  private:

    uint length;
    uint num_stored = 0;

    rep_array_unsliced<x_type> addrs;
    rep_array_unsliced<y_type> data;

  public:

    StupidLevel(uint length)
    :length(length), addrs(length), data(length)
    {
    }

    void query(rep_array_unsliced<x_type> query_addr, rep_array_unsliced<y_type> query_result, rep_array_unsliced<int> found)
    {
        uint length_for_query = max(1U, num_stored);
        const uint blocks_per_elem = blocks_per_packed_xy();
        const uint elem_stride = packed_xy_stride();
        const uint y_bytes = get_y_type_bytes();
        // Each circuit element spans blocks_per_elem blocks:
        //   layout: x_query(4B) | x(4B) | y(Y_TYPE_BYTES) | padding to block boundary
        rep_array_unsliced<block> circuit_input(length_for_query * blocks_per_elem);
        rep_array_unsliced<x_type> remove_duplicate_addr_mask(length_for_query);
        for (uint i = 0; i < length_for_query; i++) {
            circuit_input.copy_bytes_from(query_addr, sizeof(x_type), 0, i * elem_stride);
        }
        for (uint i = 0; i < num_stored; i++) {
            circuit_input.copy_bytes_from(addrs, sizeof(x_type), i * sizeof(x_type), i * elem_stride + sizeof(x_type));
            circuit_input.copy_bytes_from(data, y_bytes, i * sizeof(y_type), i * elem_stride + 2 * sizeof(x_type));
        }
        // circuit input per element:  x_query | x | y | padding
        // circuit output per element: x_mask | y_if_found | found(1b) | padding
        rep_array_unsliced<block> circuit_output(length_for_query * blocks_per_elem);
        xy_if_xs_equal_circuit->compute(circuit_output, circuit_input, length_for_query, thread_unsafe::rep_exec);

        for (uint i = 0; i < length_for_query; i++) {
            remove_duplicate_addr_mask.copy_bytes_from(circuit_output, sizeof(x_type), i * elem_stride, i * sizeof(x_type));
        }

        addrs.window(0, length_for_query).xor_with(remove_duplicate_addr_mask);

        // XOR all multi-block elements to accumulate results
        rep_array_unsliced<block> accum(blocks_per_elem);
        for (uint i = 0; i < length_for_query; i++) {
            rep_array_unsliced<block> elem = circuit_output.window(i * blocks_per_elem, blocks_per_elem);
            accum.xor_with(elem);
        }
        // Output layout: x_mask(4B) | y(Y_TYPE_BYTES) | found(1b) | padding
        query_result.copy_bytes_from(accum, y_bytes, sizeof(x_type));
        found.copy_bytes_from(accum, 1, sizeof(x_type) + y_bytes);
        accum.destroy();

        circuit_input.destroy();
        circuit_output.destroy();
        remove_duplicate_addr_mask.destroy();
    }

    void extract(rep_array_unsliced<x_type> xs, rep_array_unsliced<y_type> ys) {
        assert(num_stored == length);
        assert(xs.length_Ts() == length);
        assert(ys.length_Ts() == length);
        xs.copy(addrs);
        ys.copy(data);
    }

    void write(rep_array_unsliced<x_type> &write_addrs, rep_array_unsliced<y_type> &write_data)
    {
        assert(write_addrs.length_Ts() == write_data.length_Ts());
        assert("The stupid level is full! You should not be writing to it and have a bug elsewhere" &&
               num_stored < length);
        addrs.copy_Ts_from(write_addrs, write_addrs.length_Ts(), 0, num_stored);
        data.copy_Ts_from(write_data, write_data.length_Ts(), 0, num_stored);
        num_stored += write_addrs.length_Ts();
    }

    void skip(uint num_to_skip) {
        // this works OK with extract() because zero indices get extracted as dummies
        num_stored += num_to_skip;
    }

    bool is_writeable()
    {
        return num_stored < length;
    }

    void clear()
    {
        num_stored = 0;
    }
};
} // namespace emp

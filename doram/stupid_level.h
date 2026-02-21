#pragma once

#include "globals.h"
#include "rep_array_unsliced.h"
#include "bristol_fashion_array.h"
#include "runtime_width.h"
#include "utils.h"

#include "emp-tool/emp-tool.h"

namespace emp
{

class RuntimeStupidLevel
{
  private:
    uint length;
    uint num_stored = 0;
    BristolFashion_array* xy_if_xs_equal_circuit_local;
    runtime_width::RuntimeWidthSpec width_spec;

    rep_array_unsliced<x_type> addrs;
    rep_array_unsliced<block> data_blocks;

  public:
    RuntimeStupidLevel(
        uint length,
        BristolFashion_array* xy_if_xs_equal_circuit_local,
        const runtime_width::RuntimeWidthSpec& width_spec)
        : length(length),
          xy_if_xs_equal_circuit_local(xy_if_xs_equal_circuit_local),
          width_spec(width_spec),
          addrs(length),
          data_blocks(static_cast<uint64_t>(length) * width_spec.y_blocks) {
        assert(this->xy_if_xs_equal_circuit_local != nullptr);
        assert(this->width_spec.payload_bytes > 0);
        assert(this->width_spec.y_blocks > 0);
        assert(this->width_spec.y_stride_bytes > 0);
        assert(this->width_spec.packed_elem_blocks > 0);
        assert(this->width_spec.packed_elem_bytes > 0);
    }

    void query_blocks(
        rep_array_unsliced<x_type> query_addr,
        rep_array_unsliced<block> query_result_blocks,
        rep_array_unsliced<int> found) {
        assert(query_addr.length_Ts() == 1);
        assert(query_result_blocks.length_Ts() == width_spec.y_blocks);
        assert(found.length_Ts() == 1);
        uint length_for_query = max(1U, num_stored);
        const uint y_bytes = width_spec.payload_bytes;
        const uint blocks_per_elem = width_spec.packed_elem_blocks;
        const uint elem_stride = blocks_per_elem * sizeof(block);
        rep_array_unsliced<block> circuit_input(length_for_query * blocks_per_elem);
        rep_array_unsliced<x_type> remove_duplicate_addr_mask(length_for_query);
        for (uint i = 0; i < length_for_query; i++) {
            circuit_input.copy_bytes_from(query_addr, sizeof(x_type), 0, i * elem_stride);
        }
        for (uint i = 0; i < num_stored; i++) {
            circuit_input.copy_bytes_from(
                addrs, sizeof(x_type), i * sizeof(x_type), i * elem_stride + sizeof(x_type)
            );
            circuit_input.copy_bytes_from(
                data_blocks,
                y_bytes,
                static_cast<uint64_t>(i) * width_spec.y_stride_bytes,
                i * elem_stride + 2 * sizeof(x_type)
            );
        }
        rep_array_unsliced<block> circuit_output(length_for_query * blocks_per_elem);
        xy_if_xs_equal_circuit_local->compute(
            circuit_output, circuit_input, length_for_query, thread_unsafe::rep_exec
        );

        for (uint i = 0; i < length_for_query; i++) {
            remove_duplicate_addr_mask.copy_bytes_from(
                circuit_output, sizeof(x_type), i * elem_stride, i * sizeof(x_type)
            );
        }

        addrs.window(0, length_for_query).xor_with(remove_duplicate_addr_mask);

        rep_array_unsliced<block> accum(blocks_per_elem);
        for (uint i = 0; i < length_for_query; i++) {
            rep_array_unsliced<block> elem = circuit_output.window(i * blocks_per_elem, blocks_per_elem);
            accum.xor_with(elem);
        }
        query_result_blocks.copy_bytes_from(accum, y_bytes, sizeof(x_type));
        found.copy_bytes_from(accum, 1, sizeof(x_type) + y_bytes);
        accum.destroy();

        circuit_input.destroy();
        circuit_output.destroy();
        remove_duplicate_addr_mask.destroy();
    }

    void extract_blocks(rep_array_unsliced<x_type> xs, rep_array_unsliced<block> ys_blocks) {
        assert(num_stored == length);
        assert(xs.length_Ts() == length);
        assert(ys_blocks.length_Ts() == static_cast<uint64_t>(length) * width_spec.y_blocks);
        xs.copy(addrs);
        ys_blocks.copy(data_blocks);
    }

    void write_blocks(rep_array_unsliced<x_type>& write_addrs, rep_array_unsliced<block>& write_data_blocks) {
        const uint64_t rows = write_addrs.length_Ts();
        assert(write_data_blocks.length_Ts() == rows * width_spec.y_blocks);
        assert(
            "The stupid level is full! You should not be writing to it and have a bug elsewhere" &&
            num_stored < length
        );
        addrs.copy_Ts_from(write_addrs, rows, 0, num_stored);
        for (uint64_t i = 0; i < rows; ++i) {
            data_blocks.copy_bytes_from(
                write_data_blocks,
                width_spec.y_stride_bytes,
                i * width_spec.y_stride_bytes,
                static_cast<uint64_t>(num_stored + i) * width_spec.y_stride_bytes
            );
        }
        num_stored += rows;
    }

    void skip(uint num_to_skip) {
        num_stored += num_to_skip;
    }

    bool is_writeable() {
        return num_stored < length;
    }

    void clear() {
        num_stored = 0;
    }
};

} // namespace emp

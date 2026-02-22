#pragma once

#include <vector>
#include "emp-tool/emp-tool.h"
#include "optimal_cht.h"
#include "sh_shuffle_array.h"
#include "sh_rep_array.h"
#include "bristol_fashion_array.h"
#include "runtime_width.h"
#include "runtime_y_ops.h"
#include "utils.h"
#include "globals.h"

namespace emp
{

// couldn't put this in utils because of include order
uint prf_key_size_blocks() {
    assert(prf_circuit != nullptr);
    uint bits_per_block = 8 * sizeof(block); // 128
    assert(prf_circuit->num_input % bits_per_block == 0);
    return prf_circuit->num_input/bits_per_block - 1;
}

thread_local double time_before_cht, time_in_cht, time_after_cht, time_ite;

struct OHTableParams
{
    uint num_elements = UINT_MAX;
    uint num_dummies = UINT_MAX;
    uint stash_size = UINT_MAX;
    int builder = -1;
    uint cht_log_single_col_len = UINT_MAX;
    uint key_size_blocks = UINT_MAX;

    void validate() const {
        assert(num_elements != UINT_MAX);
        assert(num_dummies != UINT_MAX);
        assert(stash_size != UINT_MAX);
        assert(builder != -1);
        assert(cht_log_single_col_len != UINT_MAX);
    }

    uint total_size() const {
        return num_elements + num_dummies;
    }

    uint cht_full_table_length() const {
        return 2 << cht_log_single_col_len;
    }
};

ostream& operator<<(ostream& stream, const OHTableParams& params) {
    stream << "OHTable {\n";
    stream << "size: " << params.total_size() << '\n';
    stream << "num_elements: " << params.num_elements << '\n';
    stream << "num_dummies: " << params.num_dummies << '\n';
    stream << "stash_size: " << params.stash_size << '\n';
    stream << "builder: " << params.builder << '\n';
    stream << "cht_log_single_col_len: " << params.cht_log_single_col_len << '\n';
    return stream << "}\n";
}

class OHTableArray
{
    // change to private after debuggine
public:
    OHTableParams params;
    rep_array_unsliced<block> key;
    rep_array_unsliced<x_type> stash_xs;
    rep_array_unsliced<block> stash_ys_blocks;
private:
    rep_array_unsliced<block> qs_builder_order;
    rep_array_unsliced<x_type> xs_builder_order;
    rep_array_unsliced<block> ys_builder_order_blocks;
    rep_array_unsliced<uint> dummy_indices;
    rep_array_unsliced<x_type> xs_receiver_order;
    rep_array_unsliced<block> ys_receiver_order_blocks;
    block* cht_2shares;
    LocalPermutation* receiver_shuffle;
    runtime_width::RuntimeWidthSpec runtime_width_spec;
    uint64_t y_bytes;
    uint64_t y_blocks;
    uint64_t y_stride_bytes;
    uint64_t stash_y_bytes;
    uint64_t stash_y_blocks;
    uint64_t stash_y_stride_bytes;

    uint query_count = 0;

    vector<bool> touched;

public:
    OHTableArray(
        const OHTableParams& params,
        rep_array_unsliced<x_type> xs,
        rep_array_unsliced<block> ys_blocks,
        rep_array_unsliced<block> key,
        const runtime_width::RuntimeWidthSpec& runtime_spec)
    :params(params),
    key(key),
    stash_xs(params.stash_size),
    stash_ys_blocks(
        static_cast<uint64_t>(params.stash_size) * (
            (
                ((runtime_spec.total_bits + 7U) / 8U) + sizeof(block) - 1
            ) / sizeof(block)
        )
    ),
    qs_builder_order(params.total_size()),
    xs_builder_order(params.total_size()),
    ys_builder_order_blocks(
        static_cast<uint64_t>(params.total_size()) * (
            (
                ((runtime_spec.total_bits + 7U) / 8U) + sizeof(block) - 1
            ) / sizeof(block)
        )
    ),
    dummy_indices(params.num_dummies),
    xs_receiver_order(params.total_size()),
    ys_receiver_order_blocks(
        static_cast<uint64_t>(params.total_size()) * (
            (
                ((runtime_spec.total_bits + 7U) / 8U) + sizeof(block) - 1
            ) / sizeof(block)
        )
    ),
    runtime_width_spec(runtime_spec),
    y_bytes((runtime_spec.total_bits + 7U) / 8U),
    y_blocks((y_bytes + sizeof(block) - 1) / sizeof(block)),
    y_stride_bytes(y_blocks * sizeof(block)),
    stash_y_bytes((runtime_spec.total_bits + 7U) / 8U),
    stash_y_blocks((stash_y_bytes + sizeof(block) - 1) / sizeof(block)),
    stash_y_stride_bytes(stash_y_blocks * sizeof(block)),
    touched(params.total_size())
    {
        assert(xs.length_Ts() == params.num_elements);
        assert(ys_blocks.length_Ts() == static_cast<uint64_t>(params.num_elements) * y_blocks);
        assert(key.length_Ts() == prf_key_size_blocks());
        assert(y_bytes > 0);
        assert(y_blocks > 0);
        assert(y_stride_bytes > 0);
        assert(stash_y_bytes == y_bytes);
        assert(stash_y_blocks == y_blocks);
        assert(stash_y_stride_bytes == y_stride_bytes);

        params.validate();

        build_blocks(xs, ys_blocks);
    }

    ~OHTableArray() {
        key.destroy();
        stash_ys_blocks.destroy();
        qs_builder_order.destroy();
        xs_builder_order.destroy();
        ys_builder_order_blocks.destroy();
        dummy_indices.destroy();
        xs_receiver_order.destroy();
        ys_receiver_order_blocks.destroy();
    }

    void build_blocks(rep_array_unsliced<x_type> xs, rep_array_unsliced<block> ys_blocks) {
        using namespace thread_unsafe;
        keep_payload_and_alibi_blocks(ys_blocks, xs.length_Ts());
        uint prf_input_size_blocks = prf_key_size_blocks() + 1;
        rep_array_unsliced<block> keys_and_inputs(prf_input_size_blocks * params.num_elements);
        for (uint i = 0; i < params.num_elements; i++) {
            keys_and_inputs.copy_Ts_from(key, prf_key_size_blocks(), 0, prf_input_size_blocks * i);
            keys_and_inputs.copy_one<x_type>(prf_input_size_blocks * (i + 1) - 1, xs, i);
        }
        auto start = clock_start();
        prf_circuit->compute_multithreaded(qs_builder_order, keys_and_inputs, params.num_elements);
        time_total_build_prf += time_from(start);
        ArrayShuffler builder_shuffler(params.total_size());
        xs_builder_order.copy_bytes_from(xs, xs.length_bytes);
        ys_builder_order_blocks.copy_bytes_from(
            ys_blocks,
            static_cast<uint64_t>(params.num_elements) * y_stride_bytes
        );
        rep_array_unsliced<uint> indices_builder_order(params.total_size());
        builder_shuffler.forward(qs_builder_order);
        builder_shuffler.forward(xs_builder_order);
        builder_shuffler.forward_rows(ys_builder_order_blocks, static_cast<uint>(y_blocks));
        builder_shuffler.indices(indices_builder_order);
        dummy_indices.copy_Ts_from(indices_builder_order, params.num_dummies, params.num_elements);

        vector<block> qs_in_clear_compacted(params.num_elements);
        {
            vector<block> qs_in_clear(params.total_size());
            qs_builder_order.reveal_to(params.builder, qs_in_clear.data());

            if (party == params.builder) {
                uint j = 0;
                for (uint i = 0; i < params.total_size(); i++) {
                    if (blocksEqual(qs_in_clear[i], zero_block)) continue;
                    qs_in_clear_compacted[j] = (qs_in_clear[i] & makeBlock(ULLONG_MAX, ULLONG_MAX << 32)) | makeBlock(0, i);
                    j++;
                }
                assert(j == params.num_elements);
            }
        }
        vector<block> builder_cht;
        vector<uint> stash_indices_builder(params.stash_size);
        if (party == params.builder) {
            LocalPermutation builder_local_perm(thread_unsafe::private_prg, params.num_elements);
            builder_local_perm.shuffle(qs_in_clear_compacted.data());
            optimalcht::build(builder_cht, params.cht_log_single_col_len, qs_in_clear_compacted, stash_indices_builder);
        }
        rep_array_unsliced<block> cht_shares(params.cht_full_table_length());
        cht_2shares = new block[params.cht_full_table_length()];
        cht_shares.input(params.builder, builder_cht.data());
        cht_shares.reshare_3to2(prev_party(params.builder), next_party(params.builder), cht_2shares);

        ArrayShuffler receiver_shuffler(params.total_size());
        xs_receiver_order.copy(xs_builder_order);
        ys_receiver_order_blocks.copy(ys_builder_order_blocks);
        receiver_shuffler.forward_known_to_p_and_next(next_party(params.builder), xs_receiver_order);
        receiver_shuffler.forward_known_to_p_and_next_rows(
            next_party(params.builder),
            ys_receiver_order_blocks,
            static_cast<uint>(y_blocks)
        );
        if (party == prev_party(params.builder)) {
            receiver_shuffle = new LocalPermutation(receiver_shuffler.prev_shared_perm);
        } else if (party == next_party(params.builder)) {
            receiver_shuffle = new LocalPermutation(receiver_shuffler.next_shared_perm);
        }

        if (party == params.builder) {
            prev_io->send_data(stash_indices_builder.data(), params.stash_size * sizeof(uint));
            next_io->send_data(stash_indices_builder.data(), params.stash_size * sizeof(uint));
        } else if (party == next_party(params.builder)) {
            prev_io->recv_data(stash_indices_builder.data(), params.stash_size * sizeof(uint));
        } else {
            next_io->recv_data(stash_indices_builder.data(), params.stash_size * sizeof(uint));
        }

        vector<uint> stash_indices_receiver(params.stash_size);
        if (party != params.builder) {
            for (uint i = 0; i < params.stash_size; i++) {
                stash_indices_receiver[i] = receiver_shuffle->evaluate_at(stash_indices_builder[i]);
            }
        }
        
        if (party == prev_party(params.builder)) {
            next_io->send_data(stash_indices_receiver.data(), params.stash_size * sizeof(uint));
        } else if (party == params.builder) {
            prev_io->recv_data(stash_indices_receiver.data(), params.stash_size * sizeof(uint));
        }

        for (uint i = 0; i < params.stash_size; i++) {
            touched[stash_indices_receiver[i]] = true;
            stash_xs.copy_one(i, xs_receiver_order, stash_indices_receiver[i]);
            stash_ys_blocks.copy_bytes_from(
                ys_receiver_order_blocks,
                stash_y_bytes,
                static_cast<uint64_t>(stash_indices_receiver[i]) * stash_y_stride_bytes,
                static_cast<uint64_t>(i) * stash_y_stride_bytes
            );
        }

        keys_and_inputs.destroy();
        cht_shares.destroy();
    }

    void query_blocks(
        rep_array_unsliced<block> q,
        rep_array_unsliced<int> use_dummy,
        rep_array_unsliced<block> y_blocks,
        rep_array_unsliced<int> found) {
        using namespace thread_unsafe;
        assert(query_count < params.num_dummies);
        assert(q.length_Ts() == 1);
        assert(y_blocks.length_Ts() == y_blocks_count_for_rows(1));
        auto start = clock_start();
        rep_array_unsliced<block> q_or_dummy(1);
        rep_array_unsliced<block> dummy(1);
        dummy.fill_random();
        auto start_ite = clock_start();
        rep_exec->if_then_else(use_dummy, dummy, q, q_or_dummy);
        time_ite = time_from(start_ite);

        block q_clear;
        // this is the correct order, other way is ~1 round slower
        q_or_dummy.reveal_to(prev_party(params.builder), &q_clear);
        q_or_dummy.reveal_to(next_party(params.builder), &q_clear);
        if (party != params.builder) {
            assert(!blocksEqual(q_clear, zero_block));
        }

        rep_array_unsliced<uint> dummy_index(1);
        dummy_index.copy_one(0, dummy_indices, query_count);
        time_before_cht = time_from(start);
        start = clock_start();
        uint index_builder_order = optimalcht::lookup_from_2shares(cht_2shares, q_clear, params.cht_log_single_col_len, 
            dummy_index, found, params.builder);
        time_in_cht = time_from(start);

        start = clock_start();
        uint index_receiver_order;
        if (party != params.builder) {
            index_receiver_order = receiver_shuffle->evaluate_at(index_builder_order);
        }

        // waiting to receive these 4 bytes costs P1 20 mics

        if (party == prev_party(params.builder)) {
            next_io->send_data(&index_receiver_order, sizeof(uint));
        } else if (party == params.builder) {
            prev_io->recv_data(&index_receiver_order, sizeof(uint));
        }

        assert(!touched[index_receiver_order]);
        touched[index_receiver_order] = true;

        y_blocks.copy_bytes_from(
            ys_receiver_order_blocks,
            y_bytes,
            static_cast<uint64_t>(index_receiver_order) * y_stride_bytes,
            0
        );
        runtime_y_ops::keep_payload_and_alibi(
            reinterpret_cast<uint8_t*>(y_blocks.mut_prev_data()),
            runtime_width_spec
        );
        runtime_y_ops::keep_payload_and_alibi(
            reinterpret_cast<uint8_t*>(y_blocks.mut_next_data()),
            runtime_width_spec
        );
        query_count++;
        time_after_cht = time_from(start);
    }

    void extract_blocks(
        rep_array_unsliced<x_type> extract_xs,
        rep_array_unsliced<block> extract_ys_blocks) {
        assert(query_count == params.num_dummies);
        assert(extract_ys_blocks.length_Ts() == y_blocks_count_for_rows(extract_xs.length_Ts()));
        uint num_extracted = 0;
        for (uint i = 0; i < params.total_size(); i++) {
            if (touched[i]) continue;
            extract_xs.copy_one(num_extracted, xs_receiver_order, i);
            extract_ys_blocks.copy_bytes_from(
                ys_receiver_order_blocks,
                y_bytes,
                static_cast<uint64_t>(i) * y_stride_bytes,
                static_cast<uint64_t>(num_extracted) * y_stride_bytes
            );
            num_extracted++;
        }
        for (uint64_t i = 0; i < num_extracted; ++i) {
            runtime_y_ops::keep_payload_and_alibi(
                reinterpret_cast<uint8_t*>(extract_ys_blocks.mut_prev_data()) + i * y_stride_bytes,
                runtime_width_spec
            );
            runtime_y_ops::keep_payload_and_alibi(
                reinterpret_cast<uint8_t*>(extract_ys_blocks.mut_next_data()) + i * y_stride_bytes,
                runtime_width_spec
            );
        }
        assert(num_extracted == params.num_elements - params.stash_size);
    }

private:
    uint64_t y_blocks_count_for_rows(uint64_t rows) const {
        return rows * y_blocks;
    }

    inline void keep_payload_and_alibi_blocks(rep_array_unsliced<block> ys_blocks, uint64_t rows) {
        for (uint64_t i = 0; i < rows; ++i) {
            runtime_y_ops::keep_payload_and_alibi(
                reinterpret_cast<uint8_t*>(ys_blocks.mut_prev_data()) + i * y_stride_bytes,
                runtime_width_spec
            );
            runtime_y_ops::keep_payload_and_alibi(
                reinterpret_cast<uint8_t*>(ys_blocks.mut_next_data()) + i * y_stride_bytes,
                runtime_width_spec
            );
        }
    }
};

}

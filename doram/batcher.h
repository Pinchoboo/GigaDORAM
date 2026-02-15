#pragma once

#include "globals.h"
#include "rep_array_unsliced.h"
#include "bristol_fashion_array.h"

namespace emp {
namespace batcher {

// bpe = blocks per element: number of T values per logical sort element
template<typename T>
void swap_pairs(BristolFashion_array* compare_and_swap_circuit, rep_array_unsliced<T> arr,
    uint64_t chunk_size, bool invert, uint bpe = 1) {
    assert(__builtin_popcount(chunk_size) == 1 && "chunk_size power of 2");
    assert(chunk_size > 1);

    uint64_t num_elements = arr.length_Ts() / bpe;
    uint64_t even_len = (num_elements / 2) * 2;

    rep_array_unsliced<T> circuit_input(even_len * bpe);
    rep_array_unsliced<T> circuit_output(even_len * bpe);

    // ceiling(num_elements / chunk_size)
    for (uint64_t chunk_start = 0; chunk_start + chunk_size/2 < num_elements; chunk_start += chunk_size) {
        uint64_t second_half_index_upper_bound = min(chunk_size/2, num_elements - chunk_start - chunk_size/2);
        for (uint64_t i = 0; i < second_half_index_upper_bound; i++) {
            circuit_input.copy_Ts_from(arr, bpe, (chunk_start + chunk_size/2 + i) * bpe,
                (chunk_start + 2 * i + 1) * bpe);
            uint64_t first_half_index;
            // interestingly, this almost worked with the if condition backwards
            if (invert) {
                first_half_index = chunk_start + chunk_size/2 - 1 - i;
            } else {
                first_half_index = chunk_start + i;
            }
            circuit_input.copy_Ts_from(arr, bpe, first_half_index * bpe, (chunk_start + 2 * i) * bpe);
        }
    }

    compare_and_swap_circuit->compute_multithreaded(circuit_output, circuit_input, num_elements/2);

    for (uint64_t chunk_start = 0; chunk_start + chunk_size/2 < num_elements; chunk_start += chunk_size) {
        uint64_t second_half_index_upper_bound = min(chunk_size/2, num_elements - chunk_start - chunk_size/2);
        for (uint64_t i = 0; i < second_half_index_upper_bound; i++) {
            arr.copy_Ts_from(circuit_output, bpe, (chunk_start + 2 * i + 1) * bpe,
                (chunk_start + chunk_size/2 + i) * bpe);
            uint64_t first_half_index;
            if (invert) {
                first_half_index = chunk_start + chunk_size/2 - 1 - i;
            } else {
                first_half_index = chunk_start + i;
            }
            arr.copy_Ts_from(circuit_output, bpe, (chunk_start + 2 * i) * bpe, first_half_index * bpe);
        }
    }

    circuit_input.destroy();
    circuit_output.destroy();
}

template<typename T>
void butterfly_head(BristolFashion_array* compare_and_swap_circuit, rep_array_unsliced<T> arr, uint64_t chunk_size, uint bpe = 1) {
    swap_pairs(compare_and_swap_circuit, arr, chunk_size, true, bpe);
}

template<typename T>
void butterfly_body(BristolFashion_array* compare_and_swap_circuit, rep_array_unsliced<T> arr, uint64_t chunk_size, uint bpe = 1) {
    if (chunk_size == 1) return;
    swap_pairs(compare_and_swap_circuit, arr, chunk_size, false, bpe);
    butterfly_body(compare_and_swap_circuit, arr, chunk_size/2, bpe);
}

template<typename T>
void butterfly(BristolFashion_array* compare_and_swap_circuit, rep_array_unsliced<T> arr,
    uint64_t chunk_size, uint bpe = 1) {
    if (chunk_size == 1) return;
    butterfly_head(compare_and_swap_circuit, arr, chunk_size, bpe);
    butterfly_body(compare_and_swap_circuit, arr, chunk_size / 2, bpe);
}

/*
Sort consecutive chunks of length chunk_size increasing, decreasing, increasing, decreasing
*/
template<typename T>
void sort_internal(BristolFashion_array* compare_and_swap_circuit,
    rep_array_unsliced<T> arr, uint64_t chunk_size, uint bpe = 1) {
    assert(__builtin_popcount(chunk_size) == 1 && "chunk_size power of 2");
    if (chunk_size == 1) {
        return;
    }
    sort_internal(compare_and_swap_circuit, arr, chunk_size/2, bpe);
    butterfly(compare_and_swap_circuit, arr, chunk_size, bpe);
}

/*
Doesn't work on 0
*/
uint64_t least_power_of_2_greater_than_or_equal_to(uint64_t len) {
    return 1ULL << (8 * sizeof(uint64_t) - __builtin_clzll(len - 1));
}

template<typename T>
void sort (BristolFashion_array* compare_and_swap_circuit, rep_array_unsliced<T> arr, uint bpe = 1) {
    auto start = clock_start();
    uint64_t num_elements = arr.length_Ts() / bpe;
    sort_internal(compare_and_swap_circuit, arr, least_power_of_2_greater_than_or_equal_to(num_elements), bpe);
    time_total_batcher += time_from(start);
}

}
}

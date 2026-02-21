#pragma once

 // typedef ull when we use it for y_type
#include "debug.h"
#include "globals.h" //player's global prgs&ios
#include "ohtable_array.h"
#include "runtime_width.h"
#include "runtime_y_ops.h"
#include "stupid_level.h"
#include "batcher.h"

#include <fstream>     //log
#include <memory>
#include <sys/types.h> // uint

//#define q_type block;

// Hack for testing: made these editable
// TODO: proper DORAM config files
uint STASH_SIZE_USE_EMPIRICAL_CHT_BOUNDS = 8;
uint STASH_SIZE_USE_PROVEN_CHT_BOUNDS = 50;

const bool use_proven_cht_bounds = false;

namespace emp
{
class DORAM
{
  public:
    // constructor inits these
    uint log_address_space_size; 
    uint num_levels; 
    uint log_amp_factor;
    uint active_payload_bits;
    BristolFashion_array* compare_swap_circuit_local;
    rep_array_unsliced<block> prf_keys;
    runtime_width::RuntimeWidthSpec runtime_width_spec;

    uint log_sls, stupid_fill_time, stash_size; // params set by testing 
    uint amp_factor;

    double d; // ratio of CHT_col_len/actual_data_len

    // prf stuff

    //* ohtable levels are 0-indexed: 
    //* first level after stupid level is 0
    //* bottom level is num_levels - 1
    vector<OHTableArray*> ohtables;
    RuntimeStupidLevel *stupid_level;

    //?all the other ohtable params: make a struct?
    vector<uint> base_b_state_vec;

    // to log
    fstream logger;

    bool had_initial_bottom_level;
public:
    DORAM(
        uint log_address_space_size,
        uint num_levels,
        uint log_amp_factor,
        BristolFashion_array* xy_if_xs_equal_circuit_local,
        BristolFashion_array* compare_swap_circuit_local,
        uint active_payload_bits,
        rep_array_unsliced<block>* ys_no_dummy_room_blocks = nullptr)
        : 
        log_address_space_size(log_address_space_size),
        num_levels(num_levels),
        log_amp_factor(log_amp_factor),
        active_payload_bits(active_payload_bits),
        compare_swap_circuit_local(compare_swap_circuit_local),
        prf_keys(num_levels * prf_key_size_blocks()),
        runtime_width_spec{}
    {
        assert(xy_if_xs_equal_circuit_local != nullptr);
        assert(this->compare_swap_circuit_local != nullptr);
        const uint32_t total_bits = runtime_width::checked_u32_add(
            this->active_payload_bits, num_levels, "DORAM total bits"
        );
        runtime_width_spec = runtime_width::make_runtime_width_spec(
            this->active_payload_bits, num_levels, total_bits
        );
        // input should be 1 shorter than address space size
        if (ys_no_dummy_room_blocks == nullptr) {
            had_initial_bottom_level = false;
        } else {
            had_initial_bottom_level = true;
            assert(
                ys_no_dummy_room_blocks->length_Ts() ==
                static_cast<uint64_t>((1U << log_address_space_size) - 1) * runtime_width_spec.y_blocks
            );
        }

        dbg("Initializing DORAM");
        init_logger();
        decide_params(); //* ALWAYS DO THIS FIRST ALL OTHER STUFF DEPENDS ON THESE INITIALIZATIONS
        assert(num_elements_at(num_levels - 1, 1) ==
               (uint)(1 << log_address_space_size) -
                   1); // sanity check on decide params. -1 becasue we don't include 0 or 2^N
        dbg("initialized DORAM params");

        //*init stupid (runtime-width path)
        stupid_level = new RuntimeStupidLevel(
            1 << log_sls, xy_if_xs_equal_circuit_local, runtime_width_spec
        );

        //*sanity checks
        assert("must have at least 1 level (+ stupid level which is automatically included) in the DORAM!" &&
               num_levels >= 1);
        assert("must reserve N through 2N-1 as dummy locations" && log_address_space_size + 1 <= sizeof(x_type) * 8);
        // assert("test crypto version first" && use_proven_cht_bounds);

        //* doram data-specific init

        ohtables.resize(num_levels);

        if (had_initial_bottom_level) {
            //* make xs for largest level
            //*we put in 2^N-1 ys list with no room for num dummies. We move them to an array with room for num dummies
            uint bottom_level_els = num_elements_at(num_levels - 1);

            rep_array_unsliced<x_type> xs_no_dummy_room(bottom_level_els);
            vector<x_type> xs_1_to_n_clear(bottom_level_els);
            for (uint i = 0; i < bottom_level_els; i++){
                xs_1_to_n_clear[i] = i+1;
            }
            assert(xs_1_to_n_clear[0] == 1);
            assert(xs_1_to_n_clear.size() == xs_no_dummy_room.length_Ts());
            xs_no_dummy_room.input_public(xs_1_to_n_clear.data());

            dbg("creating largest level and inserting stash");
            keep_payload_only_blocks(*ys_no_dummy_room_blocks, ys_no_dummy_room_blocks->length_Ts() / runtime_width_spec.y_blocks);
            new_ohtable_of_level_blocks(num_levels - 1, xs_no_dummy_room, *ys_no_dummy_room_blocks);
            insert_stash(num_levels - 1); // this will ussuly happen as part of rebuild (the only other occusion where
        } else {
            // fake having a stash
            stupid_level->skip(stash_size);
        }
                                      // we build a level) but now we have to do it manually
        //
        //* call tests
        __all_tests_have_passed = true; //__test_rebuild_w_reinsert_driver(); //(this test is outide )

        //! I didn't deallocate ohtable after rebuild testing!
        dbg("Done building DORAM");
    }

private:

    inline void keep_payload_only_blocks(rep_array_unsliced<block>& ys_blocks, uint64_t rows = 1) {
        for (uint64_t i = 0; i < rows; ++i) {
            runtime_y_ops::keep_payload_only(
                reinterpret_cast<uint8_t*>(ys_blocks.mut_prev_data()) + i * runtime_width_spec.y_stride_bytes,
                runtime_width_spec
            );
            runtime_y_ops::keep_payload_only(
                reinterpret_cast<uint8_t*>(ys_blocks.mut_next_data()) + i * runtime_width_spec.y_stride_bytes,
                runtime_width_spec
            );
        }
    }

    inline void keep_payload_and_alibi_blocks(rep_array_unsliced<block>& ys_blocks, uint64_t rows = 1) {
        for (uint64_t i = 0; i < rows; ++i) {
            runtime_y_ops::keep_payload_and_alibi(
                reinterpret_cast<uint8_t*>(ys_blocks.mut_prev_data()) + i * runtime_width_spec.y_stride_bytes,
                runtime_width_spec
            );
            runtime_y_ops::keep_payload_and_alibi(
                reinterpret_cast<uint8_t*>(ys_blocks.mut_next_data()) + i * runtime_width_spec.y_stride_bytes,
                runtime_width_spec
            );
        }
    }

    void decide_params()
    {
        d = use_proven_cht_bounds ? 2 : 1.2;

        amp_factor = 1 << log_amp_factor;

        log_sls = log_address_space_size - (num_levels - 1) * log_amp_factor;
        //! for small input sizes only concern (can't really run large inputs now)
        stash_size = use_proven_cht_bounds ? STASH_SIZE_USE_PROVEN_CHT_BOUNDS : STASH_SIZE_USE_EMPIRICAL_CHT_BOUNDS;
        assert(stash_size < (1U << log_sls) && "stash can't be smaller than stupid level");
         
        stupid_fill_time = (1U << log_sls) - stash_size;
        
        base_b_state_vec.resize(num_levels);

        if (time_total_builds.size() < num_levels) {
            time_total_builds.resize(num_levels);
        }
    }


    void init_logger()
    {
        logger.open("doram_build_log.txt", ios::out);
        if (!logger)
        {
            cerr << "failed to open DORAM logger, exiting..." << endl;
            exit(1);
        }
        logger << "keep in mind: only p1 will write to logger\n";
    }

    void logger_write(string msg)
    {
        if (party == 1)
        {
            // dbg("in logger write", msg);
            logger << msg << endl;
        }
    }

    uint get_num_alive_levels()
    {
        uint cnt = 0;
        for (uint i = 0; i < num_levels; i++)
        {
            cnt += ohtables[i] != nullptr;
        }
        return cnt;
    }

    void new_ohtable_of_level_blocks(
        uint level_num,
        rep_array_unsliced<x_type> xs,
        rep_array_unsliced<block> ys_blocks)
    {
        auto _start = clock_start();

        assert(level_num < num_levels);
        keep_payload_and_alibi_blocks(ys_blocks, xs.length_Ts());
        if (level_num == num_levels - 1) {
            base_b_state_vec[level_num] = 1;
        } else {
            base_b_state_vec[level_num] += 1;
        }

        uint state = base_b_state_vec[level_num];
        assert(state < (uint)(1 << log_amp_factor));

        dbg("building level " + to_string(level_num) + " with state " + to_string(state));
        logger_write("building level " + to_string(level_num) + " with state " + to_string(state));

        OHTableParams params;
        params.num_elements = num_elements_at(level_num);
        params.num_dummies = get_num_dummies(level_num);
        params.stash_size = stash_size;
        params.builder = 1;
        params.cht_log_single_col_len = get_log_col_len(level_num);
        rep_array_unsliced<block> key = generate_prf_key(level_num);

        OHTableArray* new_ohtable = new OHTableArray(
            params, xs, ys_blocks, key, runtime_width_spec
        );

        time_total_builds[level_num] += time_from(_start);
        ohtables[level_num] = new_ohtable;
    }

    void delete_ohtable(uint lvl)
    {
        auto _start = clock_start();
        assert(ohtables[lvl] != nullptr);
        delete ohtables[lvl];
        ohtables[lvl] = nullptr;
        if (lvl < num_levels - 1)
        {
            if (base_b_state_vec[lvl] == amp_factor - 1) {
                base_b_state_vec[lvl] = 0;
            } else {
                // do nothing because base_b_state_vec[lvl] will be incremented in new_ohtable_of_level
            }
        }
        time_total_deletes += time_from(_start);
    }

    rep_array_unsliced<block> generate_prf_key(uint level_num) {
        rep_array_unsliced<block> key(prf_key_size_blocks());
        // TODO: replace with LowMC key generation procedure
        key.fill_random();
        prf_keys.copy_Ts_from(key, prf_key_size_blocks(), 0, level_num * prf_key_size_blocks());
        return key;
    }

    
    // we err on the larget side of things, hence the +1 for possible round-downs
    uint get_log_col_len(uint level_num, uint _state = UINT_MAX)
    {
        assert(level_num < num_levels &&
               (_state == UINT_MAX ||
                _state < (uint)(1 << log_amp_factor))); // catch obvious confusion bugs, probably will help

        uint state = (_state == UINT_MAX ? (level_num == num_levels - 1 ? 1 : base_b_state_vec[level_num]) : _state);
        uint base_b_num = level_num == num_levels - 1 ? 1 : state;
        return level_num * log_amp_factor + log_sls + 31 - __builtin_clz(d * base_b_num) +
               __builtin_popcount((d * base_b_num) != 1);
    }

    uint get_num_dummies(uint level_num) // B^i * stupid_level_size (not including stash size )
    {
        assert(level_num < num_levels); // catch obvious confusion bugs, probably will help with refactoring
        return (1 << (log_amp_factor * level_num)) * (stupid_fill_time); 
    }

    uint num_elements_at(uint level_num, uint state_override = UINT_MAX) //*this depends on the state vector
    {
        assert(level_num < num_levels); 

        if (level_num == num_levels - 1)
        {
            assert(state_override == 1 || state_override == UINT_MAX);
            //-1 because 2^N and 0 are not valid addresses
            return (1 << (log_amp_factor * level_num)) * (1 << log_sls) - 1; 
        }

        uint state = (state_override == UINT_MAX ? base_b_state_vec[level_num] : state_override);
        return (1 << (log_amp_factor * level_num)) * state * (1 << log_sls);
    }

    uint total_num_els_and_dummies(uint level_num, uint state_override = UINT_MAX) {
        return num_elements_at(level_num, state_override) + get_num_dummies(level_num);
    }

  public:
    uint get_num_levels()
    {
        return num_levels;
    } // we need this for alibi testing

    bool __all_tests_have_passed = false;
    // for now do a soft initialize which calls another function, perhaps we would want to do an
    // incremental intialize
    //? note that xs and ys need be deallocated by caller? (that's the way it is, is it a good practice?, it's
    // different
    // than @ ohtable)
    ~DORAM()
    {
        auto _start = clock_start();
        dbg("destructing DORAM");
        delete stupid_level; // stupid level is allocated in constructor and never deallocated (exccept from here)

        for (uint i = 0; i < num_levels; i++)
        {
            if (ohtables[i] != nullptr)
            {
                delete_ohtable(i);
            }
        }
        logger << "destructed DORAM succesfully\n";
        logger.close();
        dbg("destructed succesfully");
        time_total_deletes += time_from(_start);
    }

    rep_array_unsliced<block> read_and_write_blocks(
        rep_array_unsliced<x_type> qry_x,
        rep_array_unsliced<block> qry_y_blocks,
        rep_array_unsliced<int> is_write)
    {
        using namespace thread_unsafe;
        assert(is_write.length_bytes > 0);
        assert(qry_y_blocks.length_Ts() == runtime_width_spec.y_blocks);
        keep_payload_only_blocks(qry_y_blocks);

        rep_array_unsliced<int> alibi_mask(num_levels);
        rep_array_unsliced<block> y_accum_blocks(runtime_width_spec.y_blocks);
        rep_array_unsliced<int> found(1);
        rep_array_unsliced<int> use_dummy(1);

        if (!stupid_level->is_writeable())
        {
            rebuild();
        }

        auto _start = clock_start();

        uint prf_input_size_blocks = (prf_key_size_blocks() + 1);
        rep_array_unsliced<block> prf_input(prf_input_size_blocks * num_levels);
        for (uint i = 0; i < num_levels; i++) {
            if (ohtables[i] != nullptr) {
                prf_input.copy_Ts_from(prf_keys, prf_key_size_blocks(), prf_key_size_blocks() * i, prf_input_size_blocks * i);
                prf_input.copy_one(prf_input_size_blocks * (i + 1) - 1, qry_x);
            }
        }

        rep_array_unsliced<block> prf_output(num_levels);
        auto query_prf_start = clock_start();
        prf_circuit->compute(prf_output, prf_input, num_levels, rep_exec);
        time_total_query_prf += time_from(query_prf_start);

        auto query_stupid_start = clock_start();
        stupid_level->query_blocks(qry_x, y_accum_blocks, found);
        time_total_query_stupid += time_from(query_stupid_start);

        keep_payload_and_alibi_blocks(y_accum_blocks);
        extract_alibi_bits_blocks(y_accum_blocks, alibi_mask);

        for (uint i = 0; i < num_levels; i++)
        {
            if (ohtables[i] != nullptr)
            {
                use_dummy.copy_one(0, alibi_mask, i);
                use_dummy.xor_with(found);

                rep_array_unsliced<block> y_returned_blocks(runtime_width_spec.y_blocks);
                rep_array_unsliced<int> found_returned(1);
                ohtables[i]->query_blocks(prf_output.window(i, 1), use_dummy, y_returned_blocks, found_returned);
                keep_payload_and_alibi_blocks(y_returned_blocks);

                y_accum_blocks.xor_with(y_returned_blocks);
                found.xor_with(found_returned);

                keep_payload_and_alibi_blocks(y_accum_blocks);
                extract_alibi_bits_blocks(y_accum_blocks, alibi_mask);
                y_returned_blocks.destroy();
                found_returned.destroy();
            }
        }

        rep_array_unsliced<block> write_y_blocks(runtime_width_spec.y_blocks);
        for (uint i = 0; i < runtime_width_spec.y_blocks; ++i) {
            rep_exec->if_then_else(
                is_write,
                qry_y_blocks.window(i, 1),
                y_accum_blocks.window(i, 1),
                write_y_blocks.window(i, 1)
            );
        }
        keep_payload_only_blocks(write_y_blocks);
        stupid_level->write_blocks(qry_x, write_y_blocks);
        write_y_blocks.destroy();

        keep_payload_only_blocks(y_accum_blocks);

        time_total_queries += time_from(_start);
        alibi_mask.destroy();
        found.destroy();
        use_dummy.destroy();
        prf_input.destroy();
        prf_output.destroy();
        return y_accum_blocks;
    }

    void rebuild()
    {
        // stupid level must be full for us to rebuild, no point checking it
        uint rebuild_to = 0;
        bool need_to_extract_from_rebuild_to = false;
        x_type tot_num_to_extract = (uint)(1 << log_sls); // we at least gotta extract from the stupid level, we see
                                                          // how much we will need to extract with this variable/

        //*look for level we rebuild to
        for (; rebuild_to < num_levels - 1; rebuild_to++)
        {
            if (base_b_state_vec[rebuild_to] < amp_factor - 1) // once found a level that isn't fully built
            {
                if (ohtables[rebuild_to] != nullptr)
                {
                    assert(base_b_state_vec[rebuild_to] > 0);
                    need_to_extract_from_rebuild_to = true;
                    assert(num_elements_at(rebuild_to) == ohtables[rebuild_to]->params.num_elements);
                    tot_num_to_extract += num_elements_at(rebuild_to);
                }

                break; // we found the level to build, we will not collect elements from greater levels
            }
            else // if the level is fully built, we will extract all of it's contents (but it is not the level we
                 // rebuild to)
            {
                tot_num_to_extract += num_elements_at(rebuild_to);
            }
        } // if all levels (but the last level) are fully built, we would have rebuild_to incrememnted to
          // num_levels-1, and rebuild that level
        if (rebuild_to == num_levels - 1)
        {
            if (base_b_state_vec[rebuild_to] != 0) {
                assert(ohtables[rebuild_to] != nullptr);
                need_to_extract_from_rebuild_to = true;
                tot_num_to_extract += num_elements_at(rebuild_to);
            } else {
                // we actually can be in bottomless mode and have a bottom level
                assert(!had_initial_bottom_level);
                assert(ohtables[rebuild_to] == nullptr);
            }
        }

        dbg_args(base_b_state_vec, rebuild_to);
        // check that the number of elements to extract equals the number of elements
        // that will reside in the new ohtable
        assert((rebuild_to == num_levels - 1)
                || tot_num_to_extract == num_elements_at(rebuild_to, base_b_state_vec[rebuild_to] + 1));

        rep_array_unsliced<x_type> extracted_list_xs(tot_num_to_extract);
        rep_array_unsliced<block> extracted_list_ys_blocks(
            static_cast<uint64_t>(tot_num_to_extract) * runtime_width_spec.y_blocks
        );

        uint num_extracted = 0;

        //* extract_stupid

        stupid_level->extract_blocks(
            extracted_list_xs.window(0, (1 << log_sls)),
            extracted_list_ys_blocks.window(0, static_cast<uint64_t>(1 << log_sls) * runtime_width_spec.y_blocks)
        );
        num_extracted += (uint)(1 << log_sls); // there is an assert in extract checking that we onlu extract when
                                               // we extract sls many els
        // dbg("num extracted from stupid", num_extracted);
        for (uint i = 0; i < (rebuild_to + need_to_extract_from_rebuild_to); i++)
        {
            ohtables[i]->extract_blocks(
                extracted_list_xs.window(num_extracted, num_elements_at(i)),
                extracted_list_ys_blocks.window(
                    static_cast<uint64_t>(num_extracted) * runtime_width_spec.y_blocks,
                    static_cast<uint64_t>(num_elements_at(i)) * runtime_width_spec.y_blocks
                )
            );

            num_extracted += num_elements_at(i);
        }
        assert(num_extracted == tot_num_to_extract);
        keep_payload_and_alibi_blocks(extracted_list_ys_blocks, tot_num_to_extract);

        //*we now have all the xs and the ys into a list, let's re-number the dummies

        // dbg("made it right before relabel dummies for rebuild_to ", rebuild_to, " and will relabel ", num_extracted,
        //     " elements");

        if (rebuild_to == num_levels - 1) //*it is different to build to largest level
        {
            // dbg("building a largest level");
            if (had_initial_bottom_level) {
                ArrayShuffler pre_cleanse_shuffler(num_extracted);
                pre_cleanse_shuffler.forward(extracted_list_xs);
                pre_cleanse_shuffler.forward_rows(extracted_list_ys_blocks, static_cast<uint>(runtime_width_spec.y_blocks));
            }

            rep_array_unsliced<x_type> cleansed_for_bottom_level_xs(num_elements_at(num_levels - 1));
            rep_array_unsliced<block> cleansed_for_bottom_level_ys_blocks(
                static_cast<uint64_t>(num_elements_at(num_levels - 1)) * runtime_width_spec.y_blocks
            );

            cleanse_bottom_level_blocks(extracted_list_xs, extracted_list_ys_blocks, 
                cleansed_for_bottom_level_xs, cleansed_for_bottom_level_ys_blocks, 
                log_address_space_size);

            for (uint i = 0; i < num_levels; i++)
            {
                if (i == num_levels - 1 && !had_initial_bottom_level && base_b_state_vec[i] == 0) break;
                delete_ohtable(i);
            }
            new_ohtable_of_level_blocks(num_levels - 1, cleansed_for_bottom_level_xs, cleansed_for_bottom_level_ys_blocks);

            cleansed_for_bottom_level_xs.destroy();
            cleansed_for_bottom_level_ys_blocks.destroy();
        }
        else
        {
            assert(num_extracted == num_elements_at(rebuild_to, base_b_state_vec[rebuild_to] + 1));
            relabel_dummies(extracted_list_xs, log_address_space_size);
            // I suspect that we will find a problem here -- our DORAM may contain duplicates?
            for (uint i = 0; i < rebuild_to; i++)
            {
                delete_ohtable(i);
            }
            if (need_to_extract_from_rebuild_to)
            {
                delete_ohtable(rebuild_to);
            }
            new_ohtable_of_level_blocks(rebuild_to, extracted_list_xs, extracted_list_ys_blocks);
        }

        stupid_level->clear();

        // dbg("rebuild succesfull, reinserting stash...");

        insert_stash(rebuild_to);
        extracted_list_xs.destroy();
        extracted_list_ys_blocks.destroy();
    }

    void cleanse_bottom_level_blocks(rep_array_unsliced<x_type> extracted_list_xs, 
        rep_array_unsliced<block> extracted_list_ys_blocks, 
        rep_array_unsliced<x_type> cleansed_for_bottom_level_xs,
        rep_array_unsliced<block> cleansed_for_bottom_level_ys_blocks,
        uint log_N
        )
    {
        keep_payload_and_alibi_blocks(extracted_list_ys_blocks, extracted_list_xs.length_Ts());
        uint num_extracted = extracted_list_xs.length_Ts();
        rep_array_unsliced<block> cleanse_bottom_level_circuit_input(num_extracted);
        rep_array_unsliced<block> cleanse_bottom_level_circuit_output(num_extracted);
        for (uint i = 0; i < num_extracted; i++) {
            cleanse_bottom_level_circuit_input.copy_one(i, extracted_list_xs, i);
        }
        dummy_check_circuit_file[log_N]->compute_multithreaded(cleanse_bottom_level_circuit_output, cleanse_bottom_level_circuit_input, 
            num_extracted);


        if (!had_initial_bottom_level)
        {
            const uint blocks_per_elem = runtime_width_spec.packed_elem_blocks;
            const uint elem_stride = runtime_width_spec.packed_elem_blocks * sizeof(block);
            const uint y_bytes = runtime_width_spec.payload_bytes;
            // Create wider sort array: each element is BLOCKS_PER_PACKED_XY blocks
            //   layout: is_dummy(4B) | x(4B) | y(Y_TYPE_BYTES) | padding
            rep_array_unsliced<block> sort_array(num_extracted * blocks_per_elem);
            std::vector<block> zero_blocks(num_extracted * blocks_per_elem, makeBlock(0, 0));
            sort_array.input_public(zero_blocks.data());
            for (uint i = 0; i < num_extracted; i++) {
                // Copy is_dummy result from dummy_check output (first sizeof(x_type) bytes)
                sort_array.copy_bytes_from(cleanse_bottom_level_circuit_output, sizeof(x_type),
                    i * sizeof(block), i * elem_stride);
                // Copy x address
                sort_array.copy_bytes_from(extracted_list_xs, sizeof(x_type),
                    sizeof(x_type) * i, i * elem_stride + sizeof(x_type));
                // Copy y data
                sort_array.copy_bytes_from(
                    extracted_list_ys_blocks,
                    y_bytes,
                    static_cast<uint64_t>(i) * runtime_width_spec.y_stride_bytes,
                    i * elem_stride + 2 * sizeof(x_type)
                );
            }
            batcher::sort<block>(compare_swap_circuit_local, sort_array, blocks_per_elem);
            uint num_to_extract = cleansed_for_bottom_level_xs.length_Ts();
            for (uint i = 0; i < num_to_extract; i++) {
                cleansed_for_bottom_level_xs.copy_bytes_from(sort_array, sizeof(x_type),
                    i * elem_stride + sizeof(x_type), sizeof(x_type) * i);
                cleansed_for_bottom_level_ys_blocks.copy_bytes_from(
                    sort_array,
                    y_bytes,
                    i * elem_stride + 2 * sizeof(x_type),
                    static_cast<uint64_t>(i) * runtime_width_spec.y_stride_bytes
                );
            }
            keep_payload_and_alibi_blocks(cleansed_for_bottom_level_ys_blocks, num_to_extract);
            sort_array.destroy();
            relabel_dummies(cleansed_for_bottom_level_xs, log_N);
        }
        else {
            // this is a 128x inefficiency in output size
            vector<block> is_dummy(num_extracted);
            cleanse_bottom_level_circuit_output.reveal_to_all(is_dummy.data());

            uint num_real_els_found = 0;
            for (uint i = 0 ; i < num_extracted ; i++) {
                if (getLSB(is_dummy[i])) continue;
                cleansed_for_bottom_level_xs.copy_one(num_real_els_found, extracted_list_xs, i);
                cleansed_for_bottom_level_ys_blocks.copy_bytes_from(
                    extracted_list_ys_blocks,
                    runtime_width_spec.payload_bytes,
                    static_cast<uint64_t>(i) * runtime_width_spec.y_stride_bytes,
                    static_cast<uint64_t>(num_real_els_found) * runtime_width_spec.y_stride_bytes
                );
                num_real_els_found++;
            }
            // every real element in the DORAM should be collected here
            assert(num_real_els_found == cleansed_for_bottom_level_xs.length_Ts());
            keep_payload_and_alibi_blocks(cleansed_for_bottom_level_ys_blocks, num_real_els_found);
            cleanse_bottom_level_circuit_input.destroy();
            cleanse_bottom_level_circuit_output.destroy();
        }
    }

    // overwrites extracted_list_xs
    static void relabel_dummies(rep_array_unsliced<x_type> extracted_list_xs, uint log_N) {
        uint num_extracted = extracted_list_xs.length_Ts();
        rep_array_unsliced<block> relabel_dummies_circuit_input(num_extracted);
        rep_array_unsliced<block> relabel_dummies_circuit_output(num_extracted);
        rep_array_unsliced<uint> new_dummy_label(1);
        for (uint i = 0; i < num_extracted; i++) {
            relabel_dummies_circuit_input.copy_one(i, extracted_list_xs, i);
            uint new_dummy_label_clear = (1 << log_N) + i;
            new_dummy_label.input_public(&new_dummy_label_clear);
            relabel_dummies_circuit_input.copy_bytes_from(new_dummy_label, 
                sizeof(uint), 0, i * sizeof(block) + sizeof(uint));
        }
        replace_if_dummy_circuit_file[log_N]->compute_multithreaded(relabel_dummies_circuit_output, relabel_dummies_circuit_input, num_extracted);
        for (uint i = 0; i < num_extracted; i++) {
            extracted_list_xs.copy_bytes_from(relabel_dummies_circuit_output, sizeof(uint), 
                i * sizeof(block), i * sizeof(uint));
        }
        relabel_dummies_circuit_input.destroy();
        relabel_dummies_circuit_output.destroy();
        new_dummy_label.destroy();
    }

    // runs a test on rebuilding
    //? maybe I do want some params in...
    /*
        the workings of doram are so intertwined that I think pretty much no matter what I do, if I don't basically
       run DORAM, then I wouldn't be able to test it properly. For example, here, I populated a number of levels and
       fake queried them w/o reinserting and while keeping track of what I queried, and in the end (also checking
       stashes) I made sure everything made it to the rebuild_to level. Still, w/o reinsetion, the data counts would
       be wrong for the rebuild-- eh maybe not? no I think this can still work..lets try more
    */

    void insert_stash(uint level_num)
    {
        assert((level_num < num_levels) && (ohtables[level_num] != nullptr));
        OHTableArray *cur_ohtable = ohtables[level_num];

        rep_array_unsliced<x_type> stash_xs = cur_ohtable->stash_xs;
        rep_array_unsliced<block> stash_ys_blocks(stash_xs.length_Ts() * runtime_width_spec.y_blocks);
        stash_ys_blocks.copy(cur_ohtable->stash_ys_blocks);
        keep_payload_only_blocks(stash_ys_blocks, stash_xs.length_Ts());
        for (uint64_t i = 0; i < stash_xs.length_Ts(); ++i) {
            runtime_y_ops::set_alibi_bit(
                reinterpret_cast<uint8_t*>(stash_ys_blocks.mut_prev_data()) +
                    i * runtime_width_spec.y_stride_bytes,
                level_num,
                runtime_width_spec
            );
            runtime_y_ops::set_alibi_bit(
                reinterpret_cast<uint8_t*>(stash_ys_blocks.mut_next_data()) +
                    i * runtime_width_spec.y_stride_bytes,
                level_num,
                runtime_width_spec
            );
        }
        stupid_level->write_blocks(stash_xs, stash_ys_blocks);
        stash_ys_blocks.destroy();
    }

    // clears teh heirchical structure and 0's the state
    void clear_doram()
    {
        for (uint i = 0; i < num_levels; i++)
        {
            if (ohtables[i] != nullptr)
            {
                delete_ohtable(i);
            }
        }
    }

    void extract_alibi_bits_blocks(
        rep_array_unsliced<block> y_accum_blocks,
        rep_array_unsliced<int> alibi_mask) {
        for (uint i = 0; i < num_levels; i++) {
            const uint32_t bit_index = runtime_width_spec.total_bits - 1 - i;
            const bool prev_bit = runtime_y_ops::get_bit_from_y_window(
                reinterpret_cast<const uint8_t*>(y_accum_blocks.prev_data()),
                bit_index,
                runtime_width_spec
            );
            const bool next_bit = runtime_y_ops::get_bit_from_y_window(
                reinterpret_cast<const uint8_t*>(y_accum_blocks.next_data()),
                bit_index,
                runtime_width_spec
            );
            rep_array_unsliced<int> bit_share(1);
            bit_share.mut_prev_data()[0] = prev_bit ? 1 : 0;
            bit_share.mut_next_data()[0] = next_bit ? 1 : 0;
            alibi_mask.window(i, 1).xor_with(bit_share);
            bit_share.destroy();
        }
    }

    // this is the old version of this function. Because I tried to do it with no reinserting, I was short on the
    // reinserted elements for anything more than building stupid into l0, which worked
};

class RuntimeDORAMCore {
  public:
    RuntimeDORAMCore(uint active_payload_bits, uint num_levels)
        : active_payload_bits_(active_payload_bits),
          num_levels_(num_levels),
          active_payload_bytes_(active_payload_bits / 8),
          active_y_blocks_(0),
          required_total_bits_(0),
          required_bytes_(0) {
        if (active_payload_bits_ == 0 || (active_payload_bits_ % 8) != 0) {
            throw std::invalid_argument("active_payload_bits must be > 0 and divisible by 8");
        }
        required_total_bits_ = runtime_width::checked_u32_add(
            active_payload_bits_, num_levels, "RuntimeDORAM required total bits"
        );
        required_bytes_ = (required_total_bits_ + 7U) / 8U;
        active_y_blocks_ = runtime_width::y_blocks_for_total_bits(required_total_bits_);
        if (active_y_blocks_ == 0) {
            throw std::invalid_argument(
                "runtime y block count must be > 0 (payload_bits=" +
                std::to_string(active_payload_bits_) +
                ", num_levels=" + std::to_string(num_levels_) +
                ", total_bits=" + std::to_string(required_total_bits_) + ")"
            );
        }
    }

    uint active_payload_bits() const {
        return active_payload_bits_;
    }

    uint active_payload_bytes() const {
        return active_payload_bytes_;
    }

    uint active_y_blocks() const {
        return active_y_blocks_;
    }

    uint32_t required_bytes() const {
        return required_bytes_;
    }

    runtime_width::RuntimeWidthSpec runtime_width_spec() const {
        return runtime_width::make_runtime_width_spec(
            active_payload_bits_,
            num_levels_,
            required_total_bits_
        );
    }

    uint64_t expected_blocks_for_rows(uint64_t row_count, const char* ctx) const {
        return runtime_width::checked_u64_mul(
            row_count, static_cast<uint64_t>(active_y_blocks_), ctx
        );
    }

    uint64_t rows_for_block_len(uint64_t block_len, const char* ctx) const {
        if (block_len % static_cast<uint64_t>(active_y_blocks_) != 0) {
            throw std::invalid_argument(
                std::string(ctx) + ": block length does not match runtime block layout"
            );
        }
        return block_len / static_cast<uint64_t>(active_y_blocks_);
    }

    void validate_query_shape(
        const rep_array_unsliced<x_type>& qry_x,
        const rep_array_unsliced<block>& qry_y_blocks,
        const char* ctx) const {
        const uint64_t expected_blocks = expected_blocks_for_rows(
            qry_x.length_Ts(), ctx
        );
        if (qry_y_blocks.length_Ts() != expected_blocks) {
            throw std::invalid_argument("qry_y_blocks length mismatch for runtime y layout");
        }
    }

  private:
    uint active_payload_bits_;
    uint num_levels_;
    uint active_payload_bytes_;
    uint active_y_blocks_;
    uint32_t required_total_bits_;
    uint32_t required_bytes_;
};

class RuntimeDORAMEngine final {
  public:
    RuntimeDORAMEngine(
        uint log_address_space_size,
        rep_array_unsliced<block> *ys_no_dummy_room_blocks,
        uint num_levels,
        uint log_amp_factor,
        BristolFashion_array* xy_if_xs_equal_circuit_local,
        BristolFashion_array* compare_swap_circuit_local,
        const RuntimeDORAMCore& core)
        : core_(core),
          inner_(std::make_unique<DORAM>(
              log_address_space_size,
              num_levels,
              log_amp_factor,
              xy_if_xs_equal_circuit_local,
              compare_swap_circuit_local,
              core_.active_payload_bits(),
              ys_no_dummy_room_blocks)) {
    }

    rep_array_unsliced<block> read_and_write_blocks(
        rep_array_unsliced<x_type> qry_x,
        rep_array_unsliced<block> qry_y_blocks,
        rep_array_unsliced<int> is_write) {
        core_.validate_query_shape(
            qry_x, qry_y_blocks, "RuntimeDORAMEngine query block count"
        );
        return inner_->read_and_write_blocks(qry_x, qry_y_blocks, is_write);
    }

    uint get_num_levels() const {
        return inner_->get_num_levels();
    }

  private:
    RuntimeDORAMCore core_;
    std::unique_ptr<DORAM> inner_;
};

class RuntimeDORAMBackend {
  public:
    RuntimeDORAMBackend(
        uint log_address_space_size,
        rep_array_unsliced<block> *ys_no_dummy_room_blocks,
        uint num_levels,
        uint log_amp_factor,
        BristolFashion_array* xy_if_xs_equal_circuit_local,
        BristolFashion_array* compare_swap_circuit_local,
        const RuntimeDORAMCore& core)
        : inner_(std::make_unique<RuntimeDORAMEngine>(
              log_address_space_size,
              ys_no_dummy_room_blocks,
              num_levels,
              log_amp_factor,
              xy_if_xs_equal_circuit_local,
              compare_swap_circuit_local,
              core)) {
    }

    rep_array_unsliced<block> read_and_write_blocks(
        rep_array_unsliced<x_type> qry_x,
        rep_array_unsliced<block> qry_y_blocks,
        rep_array_unsliced<int> is_write) {
        return inner_->read_and_write_blocks(qry_x, qry_y_blocks, is_write);
    }

    uint get_num_levels() const {
        return inner_->get_num_levels();
    }

  private:
    std::unique_ptr<RuntimeDORAMEngine> inner_;
};

class RuntimeDORAM {
  public:
    explicit RuntimeDORAM(
        uint log_address_space_size,
        rep_array_unsliced<block> *ys_no_dummy_room_blocks,
        uint num_levels,
        uint log_amp_factor,
        BristolFashion_array* xy_if_xs_equal_circuit_local,
        BristolFashion_array* compare_swap_circuit_local,
        uint active_payload_bits)
        : core_(active_payload_bits, num_levels),
          inner_(nullptr) {
        inner_ = std::make_unique<RuntimeDORAMBackend>(
            log_address_space_size,
            ys_no_dummy_room_blocks,
            num_levels,
            log_amp_factor,
            xy_if_xs_equal_circuit_local,
            compare_swap_circuit_local,
            core_
        );
    }

    rep_array_unsliced<block> read_and_write_blocks(
        rep_array_unsliced<x_type> qry_x,
        rep_array_unsliced<block> qry_y_blocks,
        rep_array_unsliced<int> is_write) {
        core_.validate_query_shape(qry_x, qry_y_blocks, "RuntimeDORAM query block count");

        return inner_->read_and_write_blocks(qry_x, qry_y_blocks, is_write);
    }

    uint get_num_levels() {
        return inner_->get_num_levels();
    }

  private:
    RuntimeDORAMCore core_;
    std::unique_ptr<RuntimeDORAMBackend> inner_;
};

} // namespace emp

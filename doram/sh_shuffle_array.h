#include "local_permutation.h"
#include "rep_array_unsliced.h"
#include "utils.h"

namespace emp
{

class ArrayShuffler {

public:
uint len;
LocalPermutation prev_shared_perm;
LocalPermutation next_shared_perm;

ArrayShuffler (uint len)
:len(len), 
prev_shared_perm(thread_unsafe::prev_prg, len), next_shared_perm(thread_unsafe::next_prg, len)
{

}

template<typename T>
inline void forward_step(int &p, rep_array_unsliced<T> &rep_array, T* two_shares)
{
    rep_array.reshare_3to2(p, next_party(p), two_shares);
    if (party == p)
    {
        next_shared_perm.shuffle(two_shares);
    }
    else if (party == next_party(p))
    {
        prev_shared_perm.shuffle(two_shares);
    }
    rep_array.from_2shares(p, next_party(p), two_shares);
}

template<typename T>
inline void shuffle_rows(LocalPermutation& perm, T* rows, uint row_width)
{
    assert(row_width > 0);
    T* col = new T[len];
    for (uint j = 0; j < row_width; ++j) {
        for (uint i = 0; i < len; ++i) {
            col[i] = rows[i * row_width + j];
        }
        perm.shuffle(col);
        for (uint i = 0; i < len; ++i) {
            rows[i * row_width + j] = col[i];
        }
    }
    delete[] col;
}

template<typename T>
inline void inverse_shuffle_rows(LocalPermutation& perm, T* rows, uint row_width)
{
    assert(row_width > 0);
    T* col = new T[len];
    for (uint j = 0; j < row_width; ++j) {
        for (uint i = 0; i < len; ++i) {
            col[i] = rows[i * row_width + j];
        }
        perm.inverse_shuffle(col);
        for (uint i = 0; i < len; ++i) {
            rows[i * row_width + j] = col[i];
        }
    }
    delete[] col;
}

template<typename T> 
void forward(rep_array_unsliced<T> rep_array)
{
    assert(rep_array.length_Ts() == len);
    T *two_shares = new T[len];
    for (int p = 1; p <= 3; p++)
    {
        forward_step(p, rep_array, two_shares);
    }
    delete[] two_shares;
}

template<typename T>
void forward_rows(rep_array_unsliced<T> rep_array, uint row_width)
{
    assert(row_width > 0);
    assert(rep_array.length_Ts() == static_cast<uint64_t>(len) * row_width);
    T* two_shares = new T[static_cast<uint64_t>(len) * row_width];
    for (int p = 1; p <= 3; ++p)
    {
        rep_array.reshare_3to2(p, next_party(p), two_shares);
        if (party == p)
        {
            shuffle_rows(next_shared_perm, two_shares, row_width);
        }
        else if (party == next_party(p))
        {
            shuffle_rows(prev_shared_perm, two_shares, row_width);
        }
        rep_array.from_2shares(p, next_party(p), two_shares);
    }
    delete[] two_shares;
}

template<typename T>
void forward_known_to_p_and_next(int p, rep_array_unsliced<T> rep_array)
{
    assert(rep_array.length_Ts() == len);
    T *two_shares = new T[len];
    forward_step(p, rep_array, two_shares);
    delete[] two_shares;
}

template<typename T>
void forward_known_to_p_and_next_rows(int p, rep_array_unsliced<T> rep_array, uint row_width)
{
    assert(row_width > 0);
    assert(rep_array.length_Ts() == static_cast<uint64_t>(len) * row_width);
    T* two_shares = new T[static_cast<uint64_t>(len) * row_width];
    rep_array.reshare_3to2(p, next_party(p), two_shares);
    if (party == p)
    {
        shuffle_rows(next_shared_perm, two_shares, row_width);
    }
    else if (party == next_party(p))
    {
        shuffle_rows(prev_shared_perm, two_shares, row_width);
    }
    rep_array.from_2shares(p, next_party(p), two_shares);
    delete[] two_shares;
}

template<typename T> 
void inverse(rep_array_unsliced<T> rep_array)
{
    assert(rep_array.length_Ts() == len);
    T *two_shares = new T[len];
    for (int p = 3; p >= 1; p--)
    {
        rep_array.reshare_3to2(p, next_party(p), two_shares);
        if (party == p)
        {
            next_shared_perm.inverse_shuffle(two_shares);
        }
        else if (party == next_party(p))
        {
            prev_shared_perm.inverse_shuffle(two_shares);
        }
        rep_array.from_2shares(p, next_party(p), two_shares);
    }
    delete[] two_shares;
}

template<typename T>
void inverse_rows(rep_array_unsliced<T> rep_array, uint row_width)
{
    assert(row_width > 0);
    assert(rep_array.length_Ts() == static_cast<uint64_t>(len) * row_width);
    T* two_shares = new T[static_cast<uint64_t>(len) * row_width];
    for (int p = 3; p >= 1; --p)
    {
        rep_array.reshare_3to2(p, next_party(p), two_shares);
        if (party == p)
        {
            inverse_shuffle_rows(next_shared_perm, two_shares, row_width);
        }
        else if (party == next_party(p))
        {
            inverse_shuffle_rows(prev_shared_perm, two_shares, row_width);
        }
        rep_array.from_2shares(p, next_party(p), two_shares);
    }
    delete[] two_shares;
}

// T = rep_array_unsliced
void indices(rep_array_unsliced<uint> indices_)
{ 
    uint *iota = new uint[len];
    for (uint i = 0; i < len; i++) 
    {
        iota[i] = i;
    }
    indices_.input_public(iota);
    delete[] iota;
    inverse(indices_);
}

};

}

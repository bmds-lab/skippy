//
//  hash.cpp
//  indexKISS
//
//  Created by Shlomo Geva on 13/7/2023.
//

#include "hash.hpp"

/*
    MURMURHASH3
 */

uint32_t murmurHash3(uint64_t key)
{
    // hash a 64bit value to 32 bits.
    key ^= key >> 33;
    key *= 0xff51afd7ed558ccdULL;
    key ^= key >> 33;
    key *= 0xc4ceb9fe1a85ec53ULL;
    key ^= key >> 33;
    return static_cast<uint32_t>(key);
}

/*
uint32_t murmurHash3(uint64_t packedKmer) {
    // not murmurHash32, used just for testing
    uint32_t hash = static_cast<uint32_t>(packedKmer >> 32) ^ static_cast<uint32_t>(packedKmer);
    return hash;
}
*/
/*
    XOR_HASH
*/
uint32_t xorHash(uint64_t packedKmer)
{
    uint32_t hash = static_cast<uint32_t>(packedKmer >> 32) ^ static_cast<uint32_t>(packedKmer);
    return hash;
}

/*
     CANONICAL_KMER
 */
uint64_t canonical_kmer(uint64_t kmer, int KMERSIZE)
{
    return kmer ^ reverse_complement(kmer, KMERSIZE);
}

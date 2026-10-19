/*
    encode_kmer_2bit.h  (CANONICAL SHARED HEADER)

    Includes:
      - striped 16-of-21 seed encoder (pattern: 2X3X2X4X2X3)
      - reverse-complement striped encoder (same pattern)
      - Andrew Trotman encode_kmer_2bit class

    Notes:
      * STRIPED_WINDOW_SIZE is now 21 (was 20).
      * Backward-compatible wrappers are provided under the old function names
        encode_striped_3X3X4X3X3X and encode_rc_striped_3X3X4X3X3X so that
        existing code compiles without sweeping call-site renames.
*/

#pragma once

#include <cstdint>
#include <cctype>
#include <stdint.h>
#include <string>

/*
    STRIPED SEED PARAMETERS
    -----------------------
    Pattern: 2X3X2X4X2X3
    Span   : 21 bases (0..20)
    Weight : 16 bases (packed into 32 bits)
*/
static constexpr int STRIPED_WINDOW_SIZE = 21;
static constexpr int STRIPED_WEIGHT = 16;

/*
    Used offsets (0-based within the 21bp window):
      0,1, 3,4,5, 7,8, 10,11,12,13, 15,16, 18,19,20
*/
static constexpr int STRIPED_USED_INDICES[STRIPED_WEIGHT] =
{
    0, 1,        // 2
    3, 4, 5,     // 3
    7, 8,        // 2
    10,11,12,13, // 4
    15,16,       // 2
    18,19,20     // 3
};

/*
    encode_striped_2X3X2X4X2X3()
    ---------------------------
    Forward striped encoder (16 bases packed into 32 bits), rejecting non-ACGT.
    Encoding: A=0, C=1, G=2, T=3
*/
inline bool encode_striped_2X3X2X4X2X3(const char *p, uint32_t &packed)
{
    packed = 0;
    for (int i = 0; i < STRIPED_WEIGHT; ++i)
    {
        char base = (char)toupper((unsigned char)p[STRIPED_USED_INDICES[i]]);
        uint32_t code;

        switch (base)
        {
            case 'A': code = 0; break;
            case 'C': code = 1; break;
            case 'G': code = 2; break;
            case 'T': code = 3; break;
            default: return false;
        }

        packed = (packed << 2) | code;
    }
    return true;
}

/*
    encode_rc_striped_2X3X2X4X2X3()
    ------------------------------
    Reverse-complement striped encoder without constructing an RC string.

    For RC encoding, we read bases from the forward window at:
        (STRIPED_WINDOW_SIZE - 1) - STRIPED_USED_INDICES[i]
    and complement them.

    Output packing uses the same 2-bit encoding as forward:
      comp(A)=T -> 3
      comp(C)=G -> 2
      comp(G)=C -> 1
      comp(T)=A -> 0
*/
inline bool encode_rc_striped_2X3X2X4X2X3(const char *p, uint32_t &packed)
{
    packed = 0;
    for (int i = 0; i < STRIPED_WEIGHT; ++i)
    {
        const int forward_pos = (STRIPED_WINDOW_SIZE - 1) - STRIPED_USED_INDICES[i];
        char base = (char)toupper((unsigned char)p[forward_pos]);

        uint32_t code;
        switch (base)
        {
            case 'A': code = 3; break; // comp(A)=T
            case 'C': code = 2; break; // comp(C)=G
            case 'G': code = 1; break; // comp(G)=C
            case 'T': code = 0; break; // comp(T)=A
            default: return false;
        }

        packed = (packed << 2) | code;
    }
    return true;
}

/*
    Backward-compatible wrappers
    ----------------------------
    Keep the old function names so existing search/index code compiles.
    (Yes, the names no longer match the pattern; that is intentional to avoid
    a sweeping refactor during this experiment phase.)
*/
inline bool encode_striped_3X3X4X3X3X(const char *p, uint32_t &packed)
{
    return encode_striped_2X3X2X4X2X3(p, packed);
}

inline bool encode_rc_striped_3X3X4X3X3X(const char *p, uint32_t &packed)
{
    return encode_rc_striped_2X3X2X4X2X3(p, packed);
}

/*
    CLASS ENCODE_KMER_2BIT
    ----------------------
    Construct and manipulate kmers encoded in 2 bits per base (A=00, C=01, G=10, T=11)
*/
class encode_kmer_2bit
{
private:
    /*
        Fast lookup for encoding and decoding kmers
    */
    static uint64_t kmer_encoding_table[256];

public:
    /*
        ENCODE_KMER_2BIT::ENCODE_KMER_2BIT()
        ------------------------------------
    */
    encode_kmer_2bit()
    {
        kmer_encoding_table[(size_t)'A'] = 0;
        kmer_encoding_table[(size_t)'C'] = 1;
        kmer_encoding_table[(size_t)'G'] = 2;
        kmer_encoding_table[(size_t)'T'] = 3;

        kmer_encoding_table[(size_t)'a'] = 0;
        kmer_encoding_table[(size_t)'c'] = 1;
        kmer_encoding_table[(size_t)'g'] = 2;
        kmer_encoding_table[(size_t)'t'] = 3;

        kmer_encoding_table[0] = 'A';
        kmer_encoding_table[1] = 'C';
        kmer_encoding_table[2] = 'G';
        kmer_encoding_table[3] = 'T';
    }

    /*
        ENCODE_KMER_2BIT::PACK_1MER()
        -----------------------------
    */
    static uint64_t pack_1mer(char base)
    {
        return kmer_encoding_table[(size_t)base];
    }

    /*
        ENCODE_KMER_2BIT::PACK_20MER()
        ------------------------------
    */
    static uint64_t pack_20mer(const char *sequence)
    {
        uint64_t packed = 0;
        for (size_t pos = 0; pos < 20; pos++)
            packed = (packed << 2) | kmer_encoding_table[(size_t)sequence[pos]];
        return packed;
    }

    /*
        ENCODE_KMER_2BIT::PACK_32MER()
        ------------------------------
    */
    static uint64_t pack_32mer(const char *sequence)
    {
        uint64_t packed = 0;
        for (size_t pos = 0; pos < 32; pos++)
            packed = (packed << 2) | kmer_encoding_table[(size_t)sequence[pos]];
        return packed;
    }

    /*
        ENCODE_KMER_2BIT::UNPACK_20MER()
        --------------------------------
    */
    static void unpack_20mer(char *into, uint64_t packed_sequence)
    {
        for (int32_t pos = 19; pos >= 0; pos--)
            *into++ = (char)kmer_encoding_table[(packed_sequence >> (pos * 2)) & 3];
        *into = '\0';
    }

    /*
        ENCODE_KMER_2BIT::UNPACK_20MER()
        --------------------------------
    */
    static std::string unpack_20mer(uint64_t packed_sequence)
    {
        std::string sequence;
        for (int32_t pos = 19; pos >= 0; pos--)
            sequence += (char)kmer_encoding_table[(packed_sequence >> (pos * 2)) & 3];
        return sequence;
    }

    /*
        ENCODE_KMER_2BIT::REVERSE_COMPLEMENT()
        --------------------------------------
    */
    static uint64_t reverse_complement(uint64_t kmer, uint64_t bases = 20)
    {
        uint64_t result = 0;
        uint64_t complement = ~kmer;

        for (uint64_t base = 0; base < bases; base++, complement >>= 2)
            result = (result << 2) | (complement & 0x03);

        return result;
    }

    /*
        ENCODE_KMER_2BIT::REVERSE_COMPLEMENT_32MER()
        --------------------------------------------
    */
    static uint64_t reverse_complement_32mer(uint64_t kmer)
    {
        uint64_t result = 0;
        uint64_t complement = ~kmer;

        for (uint64_t base = 0; base < 32; base++, complement >>= 2)
            result = (result << 2) | (complement & 0x03);

        return result;
    }

    /*
        ENCODE_KMER_2BIT::CANONICAL_32MER()
        -----------------------------------
    */
    uint64_t canonical_32mer(uint64_t kmer)
    {
        return kmer ^ reverse_complement(kmer);
    }
};


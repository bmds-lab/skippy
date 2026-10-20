/* Encode striped kmers
 *
 *   Encode 3X3X4X3X3X 16mer to deal with high error rate in reads
 *           Copyright (c) 2025 Shlomo Geva
 *
 */

#ifndef ENCODE_STRIPED_KMER_HPP
#define ENCODE_STRIPED_KMER_HPP

#include <cctype>
#include <cstdint>

// Inline version of 3X3X4X3X3X kmer encoder (16 bases packed into 32 bits)
inline bool encode_striped_3X3X4X3X3X(const char* p, uint32_t& packed)
{
    static constexpr int used_indices[16] = {
        0,  1,  2,      // 3
        4,  5,  6,      // 3
        8,  9,  10, 11, // 4
        13, 14, 15,     // 3
        17, 18, 19      // 3
    };

    packed = 0;
    for (int i = 0; i < 16; ++i) {
        char base = toupper(p[used_indices[i]]);
        uint32_t code;
        switch (base) {
        case 'A':
            code = 0;
            break;
        case 'C':
            code = 1;
            break;
        case 'G':
            code = 2;
            break;
        case 'T':
            code = 3;
            break;
        default:
            return false; // reject kmers with ambiguous bases
        }
        packed = (packed << 2) | code;
    }
    return true;
}

// Inline version of 3X3X4X3X3X reverse-complement kmer encoder
// Applies the same stripe pattern to the RC of the 20bp window, without
// constructing the RC string.
inline bool encode_rc_striped_3X3X4X3X3X(const char* p, uint32_t& packed)
{
    static constexpr int used_indices[16] = {
        0,  1,  2,      // 3
        4,  5,  6,      // 3
        8,  9,  10, 11, // 4
        13, 14, 15,     // 3
        17, 18, 19      // 3
    };

    packed = 0;
    for (int i = 0; i < 16; ++i) {
        // In RC space we read position used_indices[i].
        // That corresponds to forward position (19 - used_indices[i]).
        char base = toupper(p[19 - used_indices[i]]);

        uint32_t code;
        switch (base) {
        case 'A':
            code = 3;
            break; // comp(A)=T -> 3
        case 'C':
            code = 2;
            break; // comp(C)=G -> 2
        case 'G':
            code = 1;
            break; // comp(G)=C -> 1
        case 'T':
            code = 0;
            break; // comp(T)=A -> 0
        default:
            return false;
        }

        packed = (packed << 2) | code;
    }
    return true;
}

#endif // ENCODE_STRIPED_KMER_HPP

/*
        ENCODE_KMER_2BIT.H
        ------------------
        Copyright (c) 2023 Andrew Trotman
*/
#pragma once
/*!
        @file
        @brief Construct and manipulate kmers encoded in 2 bits per base (A=00,
   C=01, G=10, T=11)
        @author Andrew Trotman
        @copyright 2023 Andrew Trotman
*/

#include <stdint.h>

#include <string>

/*
        CLASS ENCODE_KMER_2BIT
        ----------------------
*/
/*!
        @brief Construct and manipulate kmers encoded in 2 bits per base (A=00,
   C=01, G=10, T=11)
*/
class encode_kmer_2bit {
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
    /*!
            @brief Constructor
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
    /*!
            @brief Return the encoding of a single base.
            @param base [in] The DNA base to pack.
            @returns The packed base.
    */
    static uint64_t pack_1mer(char base)
    {
        return kmer_encoding_table[(size_t)base];
    }

    /*
            ENCODE_KMER_2BIT::PACK_20MER()
            ------------------------------
    */
    /*!
            @brief Pack a 20-mer DNA sequence (consisting of just ACGT characters)
       int a 64-bit integer
            @param sequence [in] The DNA sequence to pack.
            @returns The 20-mer packed into 64-bits
    */
    static uint64_t pack_20mer(const char* sequence)
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
    /*!
            @brief Pack a 32-mer DNA sequence (consisting of just ACGT characters)
       int a 64-bit integer
            @param sequence [in] The DNA sequence to pack.
            @returns The 32-mer packed into 64-bits
    */
    static uint64_t pack_32mer(const char* sequence)
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
    /*!
            @brief Turn a packed 20-mer into an ASCII representation (ACGT).
            @param into [out] The unpaced 20-mer (a sequence of ACTG characters,
       '\0' terminated).
            @param packed_sequence [in] The packed 20-mer.
    */
    static void unpack_20mer(char* into, uint64_t packed_sequence)
    {
        for (int32_t pos = 19; pos >= 0; pos--)
            *into++ = kmer_encoding_table[(packed_sequence >> (pos * 2)) & 3];

        *into = '\0';
    }

    /*
            ENCODE_KMER_2BIT::UNPACK_20MER()
            --------------------------------
    */
    /*!
            @brief Turn a packed 20-mer into an ASCII representation (ACGT).
            @param packed_sequence [in] The packed 20-mer.
            @returns The DNA sequence as a 20-character long string.
    */
    static std::string unpack_20mer(uint64_t packed_sequence)
    {
        std::string sequence;

        for (int32_t pos = 19; pos >= 0; pos--)
            sequence += kmer_encoding_table[(packed_sequence >> (pos * 2)) & 3];

        return sequence;
    }

    /*
            ENCODE_KMER_2BIT::REVERSE_COMPLEMENT()
            --------------------------------------
    */
    /*!
            @brief Compute the reverse complement of a packed DNA sequence.  So
       "CAT" -> "ATG" (when bases == 3).  Done without unpacking first.
            @param kmer [in] The encoded 20-mer to reverse.
            @param bases [in] The length of the DNA sequence in bases
       (default=20).
            @returns The reverse complement of the sequence.
    */
    static uint64_t reverse_complement(uint64_t kmer, uint64_t bases = 20)
    {
        uint64_t result = 0;

        /*
                Compute the complement (A<->T, C<->G), which is a bit-flip because
           A=00, T=11, C=01, G=10
        */
        uint64_t complement = ~kmer;

        /*
                Reverse the order of the bases
        */
        for (uint64_t base = 0; base < bases; base++, complement >>= 2)
            result = (result << 2) | (complement & 0x03);

        return result;
    }

    /*
            ENCODE_KMER_2BIT::REVERSE_COMPLEMENT_32MER()
            -----------------------------------------
    */
    /*!
            @brief Compute the reverse complement of a packed DNA sequence.  So
       "CAT" -> "ATG" (when bases == 3).  Done without unpacking first.
            @param kmer [in] The encoded 20-mer to reverse.
            @returns The reverse complement of the sequence.
    */
    static uint64_t reverse_complement_32mer(uint64_t kmer)
    {
        uint64_t result = 0;

        /*
                Compute the complement (A<->T, C<->G), which is a bit-flip because
           A=00, T=11, C=01, G=10
        */
        uint64_t complement = ~kmer;

        /*
                Reverse the order of the bases
        */
        for (uint64_t base = 0; base < 32; base++, complement >>= 2)
            result = (result << 2) | (complement & 0x03);

        return result;
    }

    /*
             ENCODE_KMER_2BIT::CANONICAL_32MER()
             -----------------------------------
    */
    /*!
            @brief Compute a unique key for the kmer and its reverse complement
       (done without unpacking first).
            @param kmer [in] The encoded kmer to reverse.
            @returns The canonical form
    */
    uint64_t canonical_32mer(uint64_t kmer)
    {
        return kmer ^ reverse_complement(kmer);
    }
};

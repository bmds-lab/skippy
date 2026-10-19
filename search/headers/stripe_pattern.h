/*
    stripe_pattern.h
    ----------------
    Runtime stripe pattern parser and encoder for Skippy.

    Supports arbitrary stripe patterns specified as ASCII strings, e.g.:
        "2X3X2X4X2X3"                   weight=16, window=21  (default)
        "2X3X2X4X2X3X2X3X2X4X2X3"       weight=32, window=42
        "3X4XXX2X6XX3X2"                 weight=20, window=28

    Notation:
        digits  = run of care (sampled) positions
        X's     = run of don't-care (skipped) positions

    Weight <= 32 -> packed into uint32_t
    Weight <= 64 -> packed into uint64_t
    Weight >  64 -> error

    The struct is cheap to copy and intended to be computed once at startup
    then passed by const reference everywhere.

    Optimisation (over original):
        - enc_table[256]    : replaces toupper() + switch in forward encode
        - enc_rc_table[256] : replaces toupper() + switch in RC encode
        - pos_array[64]     : replaces positions[] vector indirection in hot path
        - rc_pos_array[64]  : replaces rc_positions[] vector indirection in hot path
    All four are precomputed once in parse().
    The positions / rc_positions vectors are retained for print() and any
    non-performance-critical callers.
*/

#pragma once

#include <cstdint>
#include <cctype>
#include <cstring>
#include <string>
#include <vector>
#include <iostream>

struct StripePattern
{
    // -----------------------------------------------------------------
    // Descriptive fields (set by parse, used everywhere)
    // -----------------------------------------------------------------
    std::string      pattern_str;           // e.g. "2X3X2X4X2X3"
    int              weight      = 0;       // number of care positions
    int              window      = 0;       // total window span
    bool             use64       = false;   // true if weight > 16

    // -----------------------------------------------------------------
    // Hot-path arrays  (fixed size, cache-friendly, no pointer chasing)
    // -----------------------------------------------------------------
    int     pos_array[64]    = {};          // care positions (forward)
    int     rc_pos_array[64] = {};          // care positions (reverse complement)

    // Lookup tables: index by raw ASCII byte -> 2-bit code, or 0xFF for non-ACGT.
    // enc_table    : forward   A=0 C=1 G=2 T=3
    // enc_rc_table : RC        A=3 C=2 G=1 T=0  (complement built in)
    uint8_t enc_table[256]    = {};
    uint8_t enc_rc_table[256] = {};

    // -----------------------------------------------------------------
    // Compatibility vectors  (non-hot-path: print, diagnostics, etc.)
    // -----------------------------------------------------------------
    std::vector<int> positions;             // care positions (forward)
    std::vector<int> rc_positions;          // care positions (RC)
    std::vector<int> protect_positions;     // protect positions (forward)
    std::vector<int> rc_protect_positions;  // protect positions (RC)

    // -----------------------------------------------------------------
    // parse()
    // -----------------------------------------------------------------
    // Parses the pattern string and fills all fields.
    // Returns false and prints an error on bad input.
    // -----------------------------------------------------------------
    bool parse(const std::string &pat)
    {
        pattern_str = pat;

        positions.clear();
        rc_positions.clear();
        protect_positions.clear();
        rc_protect_positions.clear();

        weight = 0;
        window = 0;

        if (pat.empty())
        {
            std::cerr << "StripePattern: empty pattern string\n";
            return false;
        }

        // Validate: only digits and 'X'/'x'
        for (char c : pat)
        {
            if (!isdigit((unsigned char)c) && c != 'X' && c != 'x')
            {
                std::cerr << "StripePattern: invalid character '" << c
                          << "' in pattern \"" << pat << "\"\n";
                return false;
            }
        }

        // Parse alternating runs of digits (care) and X's (don't-care)
        int pos = 0;
        int i   = 0;
        const int n = (int)pat.size();

        while (i < n)
        {
            if (isdigit((unsigned char)pat[i]))
            {
                int run = 0;
                while (i < n && isdigit((unsigned char)pat[i]))
                    run = run * 10 + (pat[i++] - '0');

                if (run <= 0)
                {
                    std::cerr << "StripePattern: zero-length care run in \""
                              << pat << "\"\n";
                    return false;
                }

                for (int k = 0; k < run; ++k)
                    positions.push_back(pos++);

                weight += run;
            }
            else // 'X' or 'x'
            {
                int run = 0;
                while (i < n && (pat[i] == 'X' || pat[i] == 'x'))
                {
                    ++run;
                    ++i;
                }
                pos += run;
            }
        }

        window = pos;

        // Build protect positions
        std::vector<char> isCare(window, 0);

        for (int p : positions)
            isCare[p] = 1;

        for (int p = 0; p < window; ++p)
        {
            if (!isCare[p])
                protect_positions.push_back(p);
        }

        if (weight == 0)
        {
            std::cerr << "StripePattern: no care positions in \"" << pat << "\"\n";
            return false;
        }
        if (weight > 64)
        {
            std::cerr << "StripePattern: weight " << weight
                      << " exceeds maximum of 64 in \"" << pat << "\"\n";
            return false;
        }

        use64 = (weight > 16);

        // RC positions: read care positions right-to-left across the window,
        // stored in forward order so the encoder loops i = 0..weight-1.
        rc_positions.resize(weight);
        for (int k = 0; k < weight; ++k)
            rc_positions[k] = (window - 1) - positions[k];

        // ------------------------------------------------------------------
        // Build hot-path fixed arrays from the vectors
        // ------------------------------------------------------------------
        for (int k = 0; k < weight; ++k)
        {
            pos_array[k]    = positions[k];
            rc_pos_array[k] = rc_positions[k];
        }

        // ------------------------------------------------------------------
        // Build lookup tables
        // ------------------------------------------------------------------
        // Default: every byte is invalid (0xFF)
        memset(enc_table,    0xFF, sizeof(enc_table));
        memset(enc_rc_table, 0xFF, sizeof(enc_rc_table));

        // Forward: A=0 C=1 G=2 T=3  (upper and lower case)
        enc_table[(uint8_t)'A'] = 0;  enc_table[(uint8_t)'a'] = 0;
        enc_table[(uint8_t)'C'] = 1;  enc_table[(uint8_t)'c'] = 1;
        enc_table[(uint8_t)'G'] = 2;  enc_table[(uint8_t)'g'] = 2;
        enc_table[(uint8_t)'T'] = 3;  enc_table[(uint8_t)'t'] = 3;

        // RC: complement is built into the table  A=3 C=2 G=1 T=0
        enc_rc_table[(uint8_t)'A'] = 3;  enc_rc_table[(uint8_t)'a'] = 3;
        enc_rc_table[(uint8_t)'C'] = 2;  enc_rc_table[(uint8_t)'c'] = 2;
        enc_rc_table[(uint8_t)'G'] = 1;  enc_rc_table[(uint8_t)'g'] = 1;
        enc_rc_table[(uint8_t)'T'] = 0;  enc_rc_table[(uint8_t)'t'] = 0;

        return true;
    }

    // -----------------------------------------------------------------
    // encode32()
    // -----------------------------------------------------------------
    // Encodes a forward-strand k-mer into uint32_t (weight <= 16).
    // Returns false if any care-position base is non-ACGT.
    // -----------------------------------------------------------------
    inline bool encode32(const char *p, uint32_t &packed) const
    {
        packed = 0;
        for (int i = 0; i < weight; ++i)
        {
            const uint8_t code = enc_table[(uint8_t)p[pos_array[i]]];
            if (code == 0xFF) return false;
            packed = (packed << 2) | code;
        }
        return true;
    }

    // -----------------------------------------------------------------
    // encode64()
    // -----------------------------------------------------------------
    // Encodes a forward-strand k-mer into uint64_t (weight 17..64).
    // Returns false if any care-position base is non-ACGT.
    // -----------------------------------------------------------------
    inline bool encode64(const char *p, uint64_t &packed) const
    {
        packed = 0;
        for (int i = 0; i < weight; ++i)
        {
            const uint8_t code = enc_table[(uint8_t)p[pos_array[i]]];
            if (code == 0xFF) return false;
            packed = (packed << 2) | code;
        }
        return true;
    }

    // -----------------------------------------------------------------
    // encode_rc32()
    // -----------------------------------------------------------------
    // Encodes a reverse-complement k-mer into uint32_t (weight <= 16).
    // Complement is handled by enc_rc_table - no extra step needed.
    // Returns false if any care-position base is non-ACGT.
    // -----------------------------------------------------------------
    inline bool encode_rc32(const char *p, uint32_t &packed) const
    {
        packed = 0;
        for (int i = 0; i < weight; ++i)
        {
            const uint8_t code = enc_rc_table[(uint8_t)p[rc_pos_array[i]]];
            if (code == 0xFF) return false;
            packed = (packed << 2) | code;
        }
        return true;
    }

    // -----------------------------------------------------------------
    // encode_rc64()
    // -----------------------------------------------------------------
    // Encodes a reverse-complement k-mer into uint64_t (weight 17..64).
    // Returns false if any care-position base is non-ACGT.
    // -----------------------------------------------------------------
    inline bool encode_rc64(const char *p, uint64_t &packed) const
    {
        packed = 0;
        for (int i = 0; i < weight; ++i)
        {
            const uint8_t code = enc_rc_table[(uint8_t)p[rc_pos_array[i]]];
            if (code == 0xFF) return false;
            packed = (packed << 2) | code;
        }
        return true;
    }

    // -----------------------------------------------------------------
    // print()
    // -----------------------------------------------------------------
    void print() const
    {
        std::cout << "Stripe pattern : " << pattern_str << "\n";
        std::cout << "Weight         : " << weight << " bases\n";
        std::cout << "Window         : " << window << " bases\n";
        std::cout << "Packed type    : uint" << (use64 ? "64" : "32") << "_t\n";
        std::cout << "Care positions :";
        for (int p : positions) std::cout << " " << p;
        std::cout << "\n";
        std::cout << "RC positions   :";
        for (int p : rc_positions) std::cout << " " << p;
        std::cout << "\n";
    }
};

// =============================================================================
// DEFAULT PATTERN  (backward-compatible with existing 16-mer indexes)
// =============================================================================
static constexpr const char* DEFAULT_STRIPE_PATTERN = "2X3X2X4X2X3";

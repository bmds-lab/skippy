#ifndef GENERATE_SEED_MATRIX_HPP
#define GENERATE_SEED_MATRIX_HPP

#include <string>
#include "stripe_pattern.h"
#include "MonteCarlo.hpp"

// Generate a seed budget matrix and save to disk.
// Uses fixed count error model (fast) with greedy seed placement.
// Produces two contracts: substitution and indels.
//
// pattern          : stripe pattern
// strategy         : SKIPPY or GREEDY seed placement
// errorType        : SUBSTITUTION or INDEL
// baseName         : output filename base
// trials           : MC trials per cell (default 10000)
// maxReadLength    : maximum read length (default 2000)
// maxErrorPct      : maximum error rate % (default 35)
// targetProbability: sensitivity target (default 0.999)
//
// Output files:
//   <baseName>_greedy.bin/.csv          — substitution contract
//   <baseName>_indels_greedy.bin/.csv   — indels contract
//
// Returns true on success.
bool generateSeedMatrix(
    const StripePattern& pattern,
    PlacementStrategy    strategy,
    ErrorType            errorType,
    const std::string&   baseName,
    int                  trials            = 10000,
    int                  maxReadLength     = 2000,
    int                  maxErrorPct       = 35,
    double               targetProbability = 0.999);

#endif // GENERATE_SEED_MATRIX_HPP

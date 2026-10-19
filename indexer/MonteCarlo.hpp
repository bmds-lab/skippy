#ifndef MONTE_CARLO_HPP
#define MONTE_CARLO_HPP

#define GATE_COVERAGE 33

#include <vector>
#include "stripe_pattern.h"

enum class PlacementStrategy { SKIPPY, GREEDY };
enum class ErrorType          { SUBSTITUTION, INDEL };

// Note: ErrorModel removed — fixed count only (faster, sufficient for contract)

void initialiseGreedyCache(const StripePattern& pattern, int maxReadLength);

class MonteCarlo
{
public:
    // Compute seed start positions — precompute once per cell, reuse across trials
    static void computeSeedStarts(
        const StripePattern& pattern,
        int                  readLength,
        int                  numSeeds,
        PlacementStrategy    strategy,
        std::vector<int>&    starts);

    // Estimate P(coverage >= GATE_COVERAGE) using fixed count errors.
    // Aborts early after failThreshold failures.
    static double estimateProbability(
        const StripePattern& pattern,
        int                  readLength,
        double               errorRate,
        int                  numSeeds,
        PlacementStrategy    strategy,
        ErrorType            errorType     = ErrorType::SUBSTITUTION,
        int                  trials        = 10000,
        int                  failThreshold = 10);

    // Find minimum seed count such that P >= targetProbability.
    // hiWarm: warm start upper bound (from previous cell in row, 0 = no warm start)
    // Returns maxSeeds+1 if unachievable (sentinel).
    static int findMinimumSeeds(
        const StripePattern& pattern,
        int                  readLength,
        double               errorRate,
        double               targetProbability,
        PlacementStrategy    strategy,
        ErrorType            errorType     = ErrorType::SUBSTITUTION,
        int                  trials        = 10000,
        int                  failThreshold = 10,
        int                  hiWarm        = 0);
};

#endif // MONTE_CARLO_HPP

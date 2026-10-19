#ifndef MONTE_CARLO_HPP
#define MONTE_CARLO_HPP

#include <vector>
#include "stripe_pattern.hpp"

enum class PlacementStrategy
{
    SKIPPY,
    GREEDY
};

class MonteCarlo
{
public:

    static double estimateProbability(
        const StripePattern& pattern,
        int readLength,
        double errorRate,
        int numSeeds,
        PlacementStrategy strategy,
        int trials = 50000);

    static int findMinimumSeeds(
        const StripePattern& pattern,
        int readLength,
        double errorRate,
        double targetProbability,
        PlacementStrategy strategy,
        int trials = 50000);

private:

    static bool runTrial(
        const StripePattern& pattern,
        int readLength,
        double errorRate,
        int numSeeds,
        PlacementStrategy strategy);

    static void computeSeedStarts(
        const StripePattern& pattern,
        int readLength,
        int numSeeds,
        PlacementStrategy strategy,
        std::vector<int>& starts);
};

#endif

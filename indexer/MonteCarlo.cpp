#include "MonteCarlo.hpp"

#include <vector>
#include <random>
#include <algorithm>
#include <cmath>
#include <climits>
#include <iostream>
#include <numeric>

using namespace std;

// ── Thread-local RNG ─────────────────────────────────────────────────────────
namespace {
    thread_local mt19937_64 rng(random_device{}());
}

// ── Greedy cache ──────────────────────────────────────────────────────────────
static vector<vector<int>> greedyOrderCache;

vector<int> buildGreedyOrder(
    int readLength,
    int span,
    const vector<int>& careOffsets,
    const vector<int>& protectOffsets);

void initialiseGreedyCache(const StripePattern& pattern, int maxReadLength)
{
    greedyOrderCache.clear();
    greedyOrderCache.resize(maxReadLength + 1);

    const int span = pattern.window;
    const auto& careOffsets = pattern.positions;

    vector<bool> isCare(span, false);
    for (int off : careOffsets) isCare[off] = true;

    vector<int> protectOffsets;
    for (int i = 0; i < span; i++)
        if (!isCare[i]) protectOffsets.push_back(i);

    for (int L = span; L <= maxReadLength; ++L)
        greedyOrderCache[L] = buildGreedyOrder(L, span, careOffsets, protectOffsets);
}

// ── Seed start positions ──────────────────────────────────────────────────────
void MonteCarlo::computeSeedStarts(
    const StripePattern& pattern,
    int readLength,
    int numSeeds,
    PlacementStrategy strategy,
    vector<int>& starts)
{
    starts.clear();
    const int span = pattern.window;
    if (numSeeds <= 0) return;

    if (strategy == PlacementStrategy::SKIPPY)
    {
        int totalLength = span * numSeeds;
        int buffer = 0;

        if (totalLength + 2 * buffer <= readLength)
        {
            if (numSeeds == 1)
            {
                starts.push_back((readLength - span) / 2);
                return;
            }
            int uncoveredLength = readLength - totalLength - 2 * buffer;
            int gap = uncoveredLength / (numSeeds - 1);
            for (int i = 0; i < numSeeds; i++)
                starts.push_back(buffer + i * (span + gap));
        }
        else
        {
            if (numSeeds == 1)
            {
                starts.push_back((readLength - span) / 2);
                return;
            }
            double step = (double)(readLength - 2*buffer - span) / (double)(numSeeds - 1);
            for (int i = 0; i < numSeeds; i++)
                starts.push_back((int)(buffer + i * step));
        }
        return;
    }

    // GREEDY
    const vector<int>& ordering = greedyOrderCache[readLength];
    int n = min(numSeeds, (int)ordering.size());
    starts.insert(starts.end(), ordering.begin(), ordering.begin() + n);
}

// ── Fixed count error placement ───────────────────────────────────────────────
static void placeFixedErrors(int readLength, int numErrors, vector<char>& error)
{
    error.assign(readLength, 0);
    vector<int> positions(readLength);
    iota(positions.begin(), positions.end(), 0);
    shuffle(positions.begin(), positions.end(), rng);
    for (int i = 0; i < numErrors; i++)
        error[positions[i]] = 1;
}

static int computeNumErrors(int readLength, double errorRate, ErrorType errorType)
{
    int n = (int)round((double)readLength * errorRate);
    n = max(0, min(n, readLength));
    if (errorType == ErrorType::INDEL && (n % 2) != 0)
        n = min(n + 1, readLength); // round up to even for exact 50/50 split
    return n;
}

// ── Coverage check ────────────────────────────────────────────────────────────
// checkAllPositions: if true, seed dies if ANY position in window has error
//                   (indels sure-kill model)
//                   if false, seed dies only if a CARE position has error
//                   (substitutions model)
static bool checkCoverage(
    const StripePattern& pattern,
    const vector<int>&   starts,
    const vector<char>&  error,
    int                  readLength,
    bool                 checkAllPositions = false)
{
    const auto& offsets = pattern.positions;
    const int   window  = pattern.window;
    vector<uint64_t> coverage((readLength + 63) / 64, 0);
    int covered = 0;

    for (int start : starts)
    {
        bool survives = true;
        if (checkAllPositions)
        {
            // Indels: any error in window kills seed
            for (int off = 0; off < window; off++)
                if (error[start + off]) { survives = false; break; }
        }
        else
        {
            // Substitutions: only care positions kill seed
            for (int off : offsets)
                if (error[start + off]) { survives = false; break; }
        }
        if (!survives) continue;

        for (int off : offsets)
        {
            int pos = start + off;
            uint64_t bit = 1ULL << (pos & 63);
            uint64_t& word = coverage[pos >> 6];
            if (!(word & bit)) { word |= bit; covered++; }
        }
        if (covered >= GATE_COVERAGE) return true;
    }
    return covered >= GATE_COVERAGE;
}

// ── Core trial — fixed count errors only ─────────────────────────────────────
static bool runTrial(
    const StripePattern& pattern,
    int                  readLength,
    double               errorRate,
    const vector<int>&   starts,
    ErrorType            errorType)
{
    if (readLength <= 0 || starts.empty()) return false;

    int numErrors = computeNumErrors(readLength, errorRate, errorType);

    if (errorType == ErrorType::SUBSTITUTION)
    {
        vector<char> error;
        placeFixedErrors(readLength, numErrors, error);
        return checkCoverage(pattern, starts, error, readLength);
    }
    else // INDEL — sure kill model: any indel in window kills seed.
         // Conservative but simple — equivalent to substitutions with no
         // stripe rescue. Treats each indel position as destroying that base.
    {
        vector<char> error;
        placeFixedErrors(readLength, numErrors, error);
        return checkCoverage(pattern, starts, error, readLength, true);
    }
}

// ── estimateProbability ───────────────────────────────────────────────────────
double MonteCarlo::estimateProbability(
    const StripePattern& pattern,
    int                  readLength,
    double               errorRate,
    int                  numSeeds,
    PlacementStrategy    strategy,
    ErrorType            errorType,
    int                  trials,
    int                  failThreshold)
{
    vector<int> starts;
    computeSeedStarts(pattern, readLength, numSeeds, strategy, starts);

    int successes = 0, failures = 0;
    for (int t = 0; t < trials; t++)
    {
        if (runTrial(pattern, readLength, errorRate, starts, errorType))
            successes++;
        else if (++failures > failThreshold)
            return 0.0;
    }
    return (double)successes / (double)trials;
}

// ── findMinimumSeeds ──────────────────────────────────────────────────────────
int MonteCarlo::findMinimumSeeds(
    const StripePattern& pattern,
    int                  readLength,
    double               errorRate,
    double               targetProbability,
    PlacementStrategy    strategy,
    ErrorType            errorType,
    int                  trials,
    int                  failThreshold,
    int                  hiWarm)
{
    const int maxSeeds  = readLength - pattern.window + 1;
    if (maxSeeds <= 0) return 1; // unachievable sentinel (maxSeeds+1 = 1)

    const int minViable = max(3,
        (GATE_COVERAGE + pattern.weight - 1) / pattern.weight);

    // Binary search in [minViable, maxSeeds]
    // hiWarm is a warm start hint — start search from there if provided,
    // but always fall back to full range if needed.
    int lo = minViable;
    int hi = maxSeeds;
    int answer = -1;

    // If warm start provided, try narrow range first [lo_warm, hiWarm]
    if (hiWarm > 0 && hiWarm < maxSeeds)
    {
        int loWarm = max(minViable, (int)(hiWarm * 0.9));
        int bsLo = loWarm, bsHi = hiWarm;
        while (bsLo <= bsHi)
        {
            int mid = (bsLo + bsHi) / 2;
            vector<int> starts;
            computeSeedStarts(pattern, readLength, mid, strategy, starts);
            int s = 0, f = 0; bool ef = false;
            for (int t = 0; t < trials; t++)
            {
                if (runTrial(pattern, readLength, errorRate, starts, errorType)) s++;
                else if (++f > failThreshold) { ef = true; break; }
            }
            double p = ef ? 0.0 : (double)s / trials;
            if (p >= targetProbability) { answer = mid; bsHi = mid - 1; }
            else bsLo = mid + 1;
        }
        // If warm range found answer, try to go lower
        if (answer > minViable)
        {
            lo = minViable; hi = answer - 1;
        }
        else if (answer < 0)
        {
            // Warm range failed — search full range
            lo = minViable; hi = maxSeeds;
        }
        else
            return answer; // found at minViable
    }

    // Full binary search in [lo, hi]
    while (lo <= hi)
    {
        int mid = (lo + hi) / 2;
        vector<int> starts;
        computeSeedStarts(pattern, readLength, mid, strategy, starts);
        int s = 0, f = 0; bool ef = false;
        for (int t = 0; t < trials; t++)
        {
            if (runTrial(pattern, readLength, errorRate, starts, errorType)) s++;
            else if (++f > failThreshold) { ef = true; break; }
        }
        double p = ef ? 0.0 : (double)s / trials;
        if (p >= targetProbability) { answer = mid; hi = mid - 1; }
        else lo = mid + 1;
    }

    return (answer >= 0) ? answer : maxSeeds + 1; // maxSeeds+1 = unachievable sentinel
}

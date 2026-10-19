#include "MonteCarlo.hpp"

#include <vector>
#include <random>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <cstdlib>
#include <thread>

using namespace std;

vector<int> buildGreedyOrder(
    int readLength,
    int span,
    const vector<int>& careOffsets,
    const vector<int>& protectOffsets);

namespace
{
    thread_local mt19937_64 rng(
        random_device{}());
}

void MonteCarlo::computeSeedStarts(
    const StripePattern& pattern,
    int readLength,
    int numSeeds,
    PlacementStrategy strategy,
    std::vector<int>& starts)
{
    starts.clear();

    const int span =
        pattern.span();

    if (numSeeds <= 0)
    {
        return;
    }

    if (strategy ==
        PlacementStrategy::SKIPPY)
    {
        int totalLength =
            span * numSeeds;

        int buffer = 0;

        if (totalLength + 2 * buffer
            <= readLength)
        {
            if (numSeeds == 1)
            {
                starts.push_back(
                    (readLength - span) / 2);

                return;
            }

            int uncoveredLength =
                readLength
                - totalLength
                - 2 * buffer;

            int gap =
                uncoveredLength
                /
                (numSeeds - 1);

            for (int i = 0;
                 i < numSeeds;
                 i++)
            {
                int pos =
                    buffer
                    +
                    i * (span + gap);

                starts.push_back(
                    pos);
            }
        }
        else
        {
            if (numSeeds == 1)
            {
                starts.push_back(
                    (readLength - span) / 2);

                return;
            }

            double step =
                (double)
                (readLength
                 - 2 * buffer
                 - span)
                /
                (double)
                (numSeeds - 1);

            for (int i = 0;
                 i < numSeeds;
                 i++)
            {
                int pos =
                    (int)
                    (buffer
                     +
                     i * step);

                starts.push_back(
                    pos);
            }
        }

        return;
    }

    const vector<int>& careOffsets =
        pattern.informativeOffsets();

    vector<bool> isCare(
        span,
        false);

    for (int off : careOffsets)
    {
        isCare[off] = true;
    }

    vector<int> protectOffsets;

    for (int i = 0;
         i < span;
         i++)
    {
        if (!isCare[i])
        {
            protectOffsets.push_back(
                i);
        }
    }

    vector<int> ordering =
        buildGreedyOrder(
            readLength,
            span,
            careOffsets,
            protectOffsets);

    if ((int)ordering.size()
        < numSeeds)
    {
        numSeeds =
            (int)ordering.size();
    }

    starts.insert(
        starts.end(),
        ordering.begin(),
        ordering.begin()
        + numSeeds);
}

bool MonteCarlo::runTrial(
    const StripePattern& pattern,
    int readLength,
    double errorRate,
    int numSeeds,
    PlacementStrategy strategy)
{
    if (readLength <= 0)
    {
        return false;
    }

    if (pattern.span() > readLength)
    {
        return false;
    }

    int numErrors =
        (int)llround(
            readLength *
            errorRate);

    if (numErrors < 0)
    {
        numErrors = 0;
    }

    if (numErrors > readLength)
    {
        numErrors = readLength;
    }

    vector<char> error(
        readLength,
        0);

    vector<int> positions(
        readLength);

    for (int i = 0;
         i < readLength;
         i++)
    {
        positions[i] = i;
    }

    shuffle(
        positions.begin(),
        positions.end(),
        rng);

    for (int i = 0;
         i < numErrors;
         i++)
    {
        error[
            positions[i]
        ] = 1;
    }

    vector<int> starts;

    computeSeedStarts(
        pattern,
        readLength,
        numSeeds,
        strategy,
        starts);

    vector<uint64_t> coverage(
        (readLength + 63) / 64,
        0);

    const auto& offsets =
        pattern.informativeOffsets();

    for (int start : starts)
    {
        bool survives = true;

        for (int off : offsets)
        {
            int pos =
                start + off;

            if (error[pos])
            {
                survives = false;
                break;
            }
        }

        if (!survives)
        {
            continue;
        }

        for (int off : offsets)
        {
            int pos =
                start + off;

            coverage[
                pos >> 6
            ]
            |=
            (
                1ULL
                <<
                (pos & 63)
            );
        }
    }

    int covered = 0;

    for (uint64_t word : coverage)
    {
        covered +=
            __builtin_popcountll(
                word);
    }

    return
        covered >= 33;
}

double MonteCarlo::estimateProbability(
    const StripePattern& pattern,
    int readLength,
    double errorRate,
    int numSeeds,
    PlacementStrategy strategy,
    int trials)
{
    int success = 0;

    for (int t = 0;
         t < trials;
         t++)
    {
        if (runTrial(
                pattern,
                readLength,
                errorRate,
                numSeeds,
                strategy))
        {
            success++;
        }
    }

    return
        (double)success
        /
        (double)trials;
}

int MonteCarlo::findMinimumSeeds(
    const StripePattern& pattern,
    int readLength,
    double errorRate,
    double targetProbability,
    PlacementStrategy strategy,
    int trials)
{
    int maxSeeds =
        readLength
        -
        pattern.span()
        +
        1;

    int lo =
        (33
         +
         pattern.informativeBases()
         -
         1)
        /
        pattern.informativeBases();

    int hi =
        maxSeeds;

    int answer =
        maxSeeds;

    while (lo <= hi)
    {
        int mid =
            (lo + hi) / 2;

        double p =
            estimateProbability(
                pattern,
                readLength,
                errorRate,
                mid,
                strategy,
                trials);

        if (p >= targetProbability)
        {
            answer = mid;
            hi = mid - 1;
        }
        else
        {
            lo = mid + 1;
        }
    }

    double pMax =
        estimateProbability(
            pattern,
            readLength,
            errorRate,
            maxSeeds,
            strategy,
            trials);

    if (pMax < targetProbability)
    {
        return maxSeeds;
    }

    return answer;
}

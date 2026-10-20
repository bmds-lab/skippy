#include "buildGreedyOrder.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <set>
#include <vector>

using namespace std;

struct GreedyState {
    int readLength;
    int maxStart;

    vector<int> coverCount;
    vector<int> protectCount;

    vector<bool> protectable;

    set<int> remaining;

    vector<int> ordering;
};

static GreedyState initialiseState(int readLength, int span, const vector<int>& protectOffsets)
{
    GreedyState state;

    state.readLength = readLength;
    state.maxStart = readLength - span;

    state.coverCount.assign(readLength, 0);
    state.protectCount.assign(readLength, 0);
    state.protectable.assign(readLength, false);

    for (int start = 0; start <= state.maxStart; start++) {
        for (int off : protectOffsets) {
            int pos = start + off;

            if (pos >= 0 && pos < readLength) {
                state.protectable[pos] = true;
            }
        }
    }

    for (int s = 0; s <= state.maxStart; s++) {
        state.remaining.insert(s);
    }

    return state;
}

static bool coverageIncomplete(const GreedyState& state)
{
    for (int x : state.coverCount) {
        if (x == 0) {
            return true;
        }
    }

    return false;
}

static int computeCoverageScore(const GreedyState& state, int start, const vector<int>& careOffsets)
{
    int score = 0;

    for (int off : careOffsets) {
        int pos = start + off;

        if (state.coverCount[pos] == 0) {
            score++;
        }
    }

    return score;
}

static double computeProtectionScore(const GreedyState& state, int start,
                                     const vector<int>& protectOffsets)
{
    double score = 0.0;

    for (int off : protectOffsets) {
        int pos = start + off;

        score += log((double)state.protectCount[pos] + 1.0);
    }

    return score;
}

static int chooseBestCoverageSeed(const GreedyState& state, const vector<int>& careOffsets)
{
    int bestStart = -1;
    int bestScore = -1;

    for (int start : state.remaining) {
        int score = computeCoverageScore(state, start, careOffsets);

        if (score > bestScore) {
            bestScore = score;
            bestStart = start;
        }
    }

    return bestStart;
}

static int chooseBestProtectionSeed(const GreedyState& state, const vector<int>& protectOffsets)
{
    double bestScore = numeric_limits<double>::max();

    vector<int> tiedStarts;

    for (int start : state.remaining) {
        double score = computeProtectionScore(state, start, protectOffsets);

        if (score < bestScore) {
            bestScore = score;

            tiedStarts.clear();

            tiedStarts.push_back(start);
        }
        else if (fabs(score - bestScore) < 1e-12) {
            tiedStarts.push_back(start);
        }
    }

    return tiedStarts[tiedStarts.size() / 2];
}

static void applySeed(GreedyState& state, int start, const vector<int>& careOffsets,
                      const vector<int>& protectOffsets)
{
    state.ordering.push_back(start);

    for (int off : careOffsets) {
        state.coverCount[start + off]++;
    }

    for (int off : protectOffsets) {
        state.protectCount[start + off]++;
    }

    state.remaining.erase(start);
}

vector<int> buildGreedyOrder(int readLength, int span, const vector<int>& careOffsets,
                             const vector<int>& protectOffsets)
{
    GreedyState state = initialiseState(readLength, span, protectOffsets);

    while (!state.remaining.empty()) {
        bool coveragePhase = coverageIncomplete(state);

        int start;

        if (coveragePhase) {
            start = chooseBestCoverageSeed(state, careOffsets);
        }
        else {
            start = chooseBestProtectionSeed(state, protectOffsets);
        }
        applySeed(state, start, careOffsets, protectOffsets);
    }

    return state.ordering;
}

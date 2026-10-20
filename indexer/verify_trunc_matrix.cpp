// verify_trunc_matrix.cpp
//
// Verifies that the truncated seed matrix achieves >= targetProbability
// sensitivity for each (readLength, errorRate) cell, using MC simulation
// with exactly the same seed placement as the search.
//
// Build in /data3/skippy/Ultimate_skippyIndexer/:
//   g++ -O3 -std=c++17 -g -pthread -Iheaders \
//       -o verify_trunc_matrix \
//       verify_trunc_matrix.cpp SeedMatrix.cpp MonteCarlo.cpp \
//       buildGreedyOrder.cpp StripePattern.cpp
//
// Wait - indexer uses stripe_pattern.h not StripePattern.cpp.
// Build:
//   g++ -O3 -std=c++17 -g -pthread -Iheaders \
//       -o verify_trunc_matrix \
//       verify_trunc_matrix.cpp SeedMatrix.cpp MonteCarlo.cpp \
//       buildGreedyOrder.cpp
//
// Usage:
//   ./verify_trunc_matrix <binfile> <stripe_pattern> [trials] [target_prob]
//
// Example:
//   ./verify_trunc_matrix \
//       /data3/GRCh38/GRCh38_16_4XX4XXX1XX3X4/GRCh38_greedy.bin \
//       4XX4XXX1XX3X4 10000 0.999

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "MonteCarlo.hpp"
#include "SeedMatrix.hpp"
#include "stripe_pattern.h"

using namespace std;

int main(int argc, char* argv[])
{
    if (argc < 3) {
        cerr << "Usage: " << argv[0] << " <binfile> <stripe_pattern> [trials] [target_prob]\n";
        return 1;
    }

    string binFile = argv[1];
    string patternStr = argv[2];
    int trials = (argc > 3) ? atoi(argv[3]) : 10000;
    double targetProb = (argc > 4) ? atof(argv[4]) : 0.999;

    // Load seed matrix
    SeedMatrix sm;
    if (!sm.load(binFile)) {
        cerr << "Failed to load " << binFile << "\n";
        return 1;
    }

    sm.printSummary();

    if (!sm.hasTruncMatrix()) {
        cerr << "No truncation matrix in " << binFile << " — cannot verify.\n";
        return 1;
    }

    // Parse stripe pattern
    StripePattern pattern;
    if (!pattern.parse(patternStr)) {
        cerr << "Failed to parse pattern: " << patternStr << "\n";
        return 1;
    }

    // Initialise greedy cache for all read lengths
    cout << "Initialising greedy order table...\n";
    initialiseGreedyCache(pattern, 1000);

    const auto& readLengths = sm.getReadLengths();
    const auto& errorRates = sm.getErrorRates();
    const int rows = (int)errorRates.size();
    const int cols = (int)readLengths.size();
    const int span = pattern.window;

    int totalCells = 0;
    int failCells = 0;
    int skipCells = 0;
    double sumAchieved = 0.0;
    double minAchieved = 1.0;
    int minR = 0, minC = 0;

    // Print header
    cout << fixed << setprecision(4);
    cout << "\nVerifying truncated matrix: " << patternStr << "\n";
    cout << "Trials: " << trials << "  Target: " << targetProb << "\n\n";

    cout << left << setw(12) << "ReadLen" << right << setw(10) << "ErrorRate" << right << setw(12)
         << "TruncLen" << right << setw(10) << "Seeds" << right << setw(10) << "Achieved" << right
         << setw(8) << "Pass"
         << "\n";
    cout << string(62, '-') << "\n";

    for (int r = 0; r < rows; r++) {
        for (int c = 0; c < cols; c++) {
            int l = readLengths[c];
            double e = errorRates[r];

            // Query truncated matrix — same as search does
            auto ts = sm.queryTrunc(e, l);

            int truncLen = ts.truncLength;
            int seedCount = ts.seedCount;
            int maxSeeds = truncLen - span + 1;

            // Skip unsatisfiable cells
            if (seedCount > maxSeeds || maxSeeds <= 0) {
                skipCells++;
                continue;
            }

            // Run MC with greedy placement on truncLen, seedCount seeds
            // This mirrors exactly what the search does
            double achieved =
                MonteCarlo::estimateProbability(pattern,
                                                truncLen, // use truncated length
                                                e, seedCount, PlacementStrategy::GREEDY, trials);

            totalCells++;
            sumAchieved += achieved;

            bool pass = (achieved >= targetProb);
            if (!pass) {
                failCells++;
                if (achieved < minAchieved) {
                    minAchieved = achieved;
                    minR = r;
                    minC = c;
                }
            }

            // Print every 50bp and every error rate, or failures
            if (!pass || (l % 50 == 0)) {
                cout << left << setw(12) << l << right << setw(9) << (int)round(e * 100) << "%"
                     << right << setw(12) << truncLen << right << setw(10) << seedCount << right
                     << setw(10) << achieved << right << setw(8) << (pass ? "PASS" : "FAIL")
                     << "\n";
            }
        }
    }

    cout << "\n" << string(62, '=') << "\n";
    cout << "Summary:\n";
    cout << "  Total verified cells : " << totalCells << "\n";
    cout << "  Skipped (unsat)      : " << skipCells << "\n";
    cout << "  Passed               : " << (totalCells - failCells) << "\n";
    cout << "  Failed               : " << failCells << "\n";

    if (totalCells > 0) {
        cout << "  Mean achieved        : " << sumAchieved / totalCells << "\n";
    }

    if (failCells > 0) {
        cout << "  Worst cell           : l=" << readLengths[minC]
             << " e=" << (int)round(errorRates[minR] * 100) << "%"
             << " achieved=" << minAchieved << "\n";
    }
    else {
        cout << "  All cells PASSED at target " << targetProb << "\n";
    }

    return failCells > 0 ? 1 : 0;
}

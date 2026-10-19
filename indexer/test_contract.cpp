// test_contract.cpp — standalone test harness for contract generation
//
// Usage:
//   ./test_contract [pattern] [placement] [error_type] [trials] [max_read] [max_error_pct]
//
// Defaults:
//   pattern       = 2X3X2X4X2X3
//   placement     = skippy         (or greedy)
//   error_type    = subs           (or indels)
//   trials        = 1000
//   max_read      = 500
//   max_error_pct = 25
//
// Outputs:
//   seed_matrix_subs_skippy.bin + .csv    (or appropriate suffix)
//
// Example:
//   ./test_contract 2X3X2X4X2X3 greedy subs 2000 1000 30
//   ./test_contract 2X3X2X4X2X3 skippy indels 1000 500 20

#include "generate_seed_matrix.hpp"
#include "MonteCarlo.hpp"
#include "stripe_pattern.h"

#include <iostream>
#include <cstring>
#include <cstdlib>
#include <string>

using namespace std;

int main(int argc, char* argv[])
{
    // ── Defaults ──────────────────────────────────────────────────────────
    string patternStr  = "2X3X2X4X2X3";
    string placementStr = "skippy";
    string errorTypeStr = "subs";
    int    trials       = 1000;
    int    maxRead      = 500;
    int    maxErrorPct  = 25;

    // ── Parse arguments ───────────────────────────────────────────────────
    if (argc > 1) patternStr   = argv[1];
    if (argc > 2) placementStr = argv[2];
    if (argc > 3) errorTypeStr = argv[3];
    if (argc > 4) trials       = atoi(argv[4]);
    if (argc > 5) maxRead      = atoi(argv[5]);
    if (argc > 6) maxErrorPct  = atoi(argv[6]);

    // ── Validate and map placement ────────────────────────────────────────
    PlacementStrategy strategy;
    if (placementStr == "greedy")
        strategy = PlacementStrategy::GREEDY;
    else if (placementStr == "skippy")
        strategy = PlacementStrategy::SKIPPY;
    else
    {
        cerr << "ERROR: placement must be 'skippy' or 'greedy'\n";
        return 1;
    }

    // ── Validate and map error type ───────────────────────────────────────
    ErrorType errorType;
    if (errorTypeStr == "subs" || errorTypeStr == "substitution")
        errorType = ErrorType::SUBSTITUTION;
    else if (errorTypeStr == "indels" || errorTypeStr == "indel")
        errorType = ErrorType::INDEL;
    else
    {
        cerr << "ERROR: error_type must be 'subs' or 'indels'\n";
        return 1;
    }

    // ── Parse stripe pattern ──────────────────────────────────────────────
    StripePattern pattern;
    if (!pattern.parse(patternStr))
    {
        cerr << "ERROR: could not parse pattern '" << patternStr << "'\n";
        return 1;
    }

    // ── Print configuration ───────────────────────────────────────────────
    cout << "=== Contract generation test ===\n"
         << "Pattern:      " << patternStr << "\n"
         << "  window:     " << pattern.window << "\n"
         << "  weight:     " << pattern.weight << "\n"
         << "Placement:    " << placementStr << "\n"
         << "Error type:   " << errorTypeStr << "\n"
         << "Trials:       " << trials << "\n"
         << "Max read len: " << maxRead << "bp\n"
         << "Max error:    " << maxErrorPct << "%\n"
         << "Target P:     0.999\n"
         << "Output base:  seed_matrix\n\n";

    // ── Generate contract ─────────────────────────────────────────────────
    bool ok = generateSeedMatrix(
        pattern,
        strategy,
        errorType,
        "seed_matrix",    // output base name
        trials,
        maxRead,
        maxErrorPct,
        0.999);

    if (!ok)
    {
        cerr << "Contract generation FAILED\n";
        return 1;
    }

    cout << "Done. Check seed_matrix*.csv for results.\n";
    return 0;
}

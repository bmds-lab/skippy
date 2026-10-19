// report_seed_matrix.cpp
#include "report_seed_matrix.hpp"
#include "SeedMatrix.hpp"

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <climits>
#include <vector>
#include <string>
#include <algorithm>

using namespace std;

// ---------------------------------------------------------------------------
// CSV helper
// ---------------------------------------------------------------------------
static void saveCSV(
    const string&           path,
    const vector<uint16_t>& data,
    const vector<int>&      readLengths,
    const vector<double>&   errorRates)
{
    const int rows = (int)errorRates.size();
    const int cols = (int)readLengths.size();
    ofstream f(path);
    f << "ErrorRate";
    for (int w : readLengths) f << "," << w;
    f << "\n";
    for (int r = 0; r < rows; r++)
    {
        f << errorRates[r];
        for (int c = 0; c < cols; c++)
            f << "," << data[r * cols + c];
        f << "\n";
    }
    cout << "Saved " << path << "\n";
}

// ---------------------------------------------------------------------------
// Main report function
// ---------------------------------------------------------------------------
void reportSeedMatrix(
    SeedMatrix&                      sm,
    const std::vector<uint16_t>&     matrix,
    const std::vector<int>&          readLengths,
    const std::vector<double>&       errorRates,
    const StripePattern&             pattern,
    const std::string&               baseName,
    PlacementStrategy                strategy)
{
    const int rows = (int)errorRates.size();
    const int cols = (int)readLengths.size();
    const int span = pattern.window;

    // ── Build truncation matrix ───────────────────────────────────────────
    // Single left-to-right pass per row.
    // Track dominant cell — shortest viable read with lowest budget.
    // Unachievable cells (matrix[idx] > maxS) store maxS.
    vector<uint16_t> truncLen      (rows * cols, 0);
    vector<uint16_t> truncSeeds    (rows * cols, 0);
    vector<uint16_t> truncErrorRate(rows * cols, 0);

    for (int r = 0; r < rows; r++)
    {
        int dominantLen   = -1;
        int dominantSeeds = INT_MAX;

        for (int c = 0; c < cols; c++)
        {
            int idx  = r * cols + c;
            int maxS = readLengths[c] - span + 1;
            int v    = (int)matrix[idx];

            if (v > maxS)  // unachievable
            {
                truncLen      [idx] = (uint16_t)readLengths[c];
                truncSeeds    [idx] = (uint16_t)maxS;
                truncErrorRate[idx] = (uint16_t)rows; // sentinel
                continue;
            }

            if (v < dominantSeeds)
            {
                dominantSeeds = v;
                dominantLen   = readLengths[c];
            }

            truncLen      [idx] = (uint16_t)dominantLen;
            truncSeeds    [idx] = (uint16_t)dominantSeeds;
            truncErrorRate[idx] = (uint16_t)r;
        }
    }

    // ── Store in SeedMatrix and save bin ──────────────────────────────────
    sm.setTruncLengths   (truncLen);
    sm.setTruncSeeds     (truncSeeds);
    sm.setTruncErrorRates(truncErrorRate);

    string stratStr = (strategy == PlacementStrategy::GREEDY) ? "_greedy" : "_skippy";
    string binName  = baseName + stratStr + ".bin";

    if (!sm.save(binName))
        cerr << "Failed to write truncation data to " << binName << "\n";
    else
        cout << "Updated " << binName << " with truncation matrix (version 3)\n";

    // ── Save CSVs ─────────────────────────────────────────────────────────
    saveCSV(baseName + stratStr + "_trunc_len.csv",
            truncLen,     readLengths, errorRates);
    saveCSV(baseName + stratStr + "_trunc_seeds.csv",
            truncSeeds,   readLengths, errorRates);
    saveCSV(baseName + stratStr + "_trunc_errrate.csv",
            truncErrorRate, readLengths, errorRates);

    // ── Contract report ───────────────────────────────────────────────────
    int    satCells          = 0;
    int    truncBenefitCells = 0;
    double sumSeeds          = 0.0;
    double sumSeedsSquared   = 0.0;

    for (int r = 0; r < rows; r++)
        for (int c = 0; c < cols; c++)
        {
            int idx = r * cols + c;
            if ((int)truncErrorRate[idx] == r)
            {
                satCells++;
                int n = (int)truncSeeds[idx];
                sumSeeds       += n;
                sumSeedsSquared += (double)n * n;
                if ((int)truncLen[idx] < readLengths[c]) truncBenefitCells++;
            }
        }

    double meanSeeds = satCells > 0 ? sumSeeds / satCells : 0.0;
    double rmsSeeds  = satCells > 0 ? sqrt(sumSeedsSquared / satCells) : 0.0;

    string reportPath = baseName + stratStr + "_contract.txt";
    ofstream report(reportPath);

    ios_base::fmtflags savedFlags = cout.flags();
    streamsize savedPrec = cout.precision();

    auto printReport = [&](ostream& out)
    {
        out << "======================================================\n"
            << "Skippy Contract Report: " << pattern.pattern_str << "\n"
            << "Strategy: " << stratStr.substr(1) << "\n"
            << "======================================================\n\n"
            << "Summary (satisfiable cells only):\n"
            << "  Grid size            : " << rows << " x " << cols
            << " = " << rows*cols << " cells\n"
            << fixed << setprecision(3)
            << "  Satisfiable fraction : "
            << 100.0 * satCells / (rows*cols) << "%"
            << "  (" << satCells << " / " << rows*cols << ")\n"
            << "  Mean seed count      : " << meanSeeds << "\n"
            << "  RMS seed count       : " << rmsSeeds  << "\n\n"
            << "Truncation savings:\n"
            << "  Cells benefiting     : " << truncBenefitCells
            << "  (" << fixed << setprecision(1)
            << 100.0 * truncBenefitCells / max(1, satCells) << "%)\n\n"
            << "Contract boundary (highest satisfiable error rate per read length):\n"
            << "  Read length     Max error rate  Seeds required"
            << "    Trunc length   Trunc seeds\n"
            << "  " << string(76, '-') << "\n";

        for (int c = 0; c < cols; c++)
        {
            if (cols > 20 && readLengths[c] % 50 != 0) continue;

            int cliffRow = -1;
            for (int r = rows-1; r >= 0; r--)
                if ((int)truncErrorRate[r*cols+c] == r)
                { cliffRow = r; break; }
            if (cliffRow < 0) continue;

            int idx     = cliffRow * cols + c;
            int tlen    = (int)truncLen[idx];
            int tseeds  = (int)truncSeeds[idx];
            int orig    = (int)matrix[cliffRow*cols+c];
            bool trunc  = (tlen < readLengths[c]);

            out << "  " << setw(11) << readLengths[c]
                << setw(20)
                << (to_string((int)round(errorRates[cliffRow]*100)) + "%")
                << setw(17) << orig
                << setw(15) << tlen << (trunc ? " *" : "  ")
                << setw(12) << tseeds << "\n";
        }
        out << "  (* = truncation beneficial)\n";
    };

    printReport(cout);
    printReport(report);
    cout << "\nContract report written to " << reportPath << "\n";

    cout.flags(savedFlags);
    cout.precision(savedPrec);
}

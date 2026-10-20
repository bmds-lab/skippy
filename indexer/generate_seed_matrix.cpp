#include "generate_seed_matrix.hpp"

#include <algorithm>
#include <atomic>
#include <climits>
#include <cmath>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <mutex>
#include <numeric>
#include <sstream>
#include <thread>
#include <vector>

#include "MonteCarlo.hpp"
#include "SeedMatrix.hpp"
#include "report_seed_matrix.hpp"
#include "stripe_pattern.h"

using namespace std;

// ─────────────────────────────────────────────────────────────────────────────
// Monotone non-increasing envelope (left to right).
// Longer reads need <= seeds than shorter reads.
// Unachievable cells (budget == maxSeeds+1 sentinel) are left unchanged.
// ─────────────────────────────────────────────────────────────────────────────
static void applyMonotoneEnvelope(vector<uint16_t>& row, const vector<int>& readLengths, int window)
{
    const int cols = (int)row.size();
    uint16_t minSeen = 0;
    bool found = false;

    for (int c = 0; c < cols; c++) {
        int maxS = readLengths[c] - window + 1;
        if ((int)row[c] > maxS)
            continue; // unachievable or pending — skip

        if (!found) {
            minSeen = row[c];
            found = true;
        }
        else if (row[c] < minSeen)
            minSeen = row[c];
        else
            row[c] = minSeen;
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Main entry point
// ─────────────────────────────────────────────────────────────────────────────
bool generateSeedMatrix(const StripePattern& pattern, PlacementStrategy strategy,
                        ErrorType errorType, const string& baseName, int trials, int maxReadLength,
                        int maxErrorPct, double targetProbability)
{
    const int minReadLength = 60;
    const int failThreshold = max(1, (int)((1.0 - targetProbability) * trials));

    vector<int> readLengths;
    for (int w = minReadLength; w <= maxReadLength; w++)
        readLengths.push_back(w);

    vector<double> errorRates;
    for (int e = 0; e <= maxErrorPct; e++)
        errorRates.push_back((double)e / 100.0);

    const int rows = (int)errorRates.size();
    const int cols = (int)readLengths.size();

    // Suffix: substitution keeps original name for backward compatibility
    string suffix = string(errorType == ErrorType::INDEL ? "_indels" : "") +
                    string(strategy == PlacementStrategy::GREEDY ? "_greedy" : "_skippy");

    {
        // Restore cout format after any previous fixed/setprecision calls
        ios_base::fmtflags savedFlags = cout.flags();
        streamsize savedPrec = cout.precision();

        cout << "\nGenerating " << (errorType == ErrorType::INDEL ? "INDEL" : "SUBSTITUTION")
             << " contract (" << (strategy == PlacementStrategy::GREEDY ? "GREEDY" : "SKIPPY")
             << ", fixed count)\n"
             << "  Read lengths: " << minReadLength << " to " << maxReadLength << "\n"
             << "  Error rates:  0% to " << maxErrorPct << "%\n"
             << fixed << setprecision(4) << "  Target P:     " << targetProbability << "\n";

        cout.flags(savedFlags);
        cout.precision(savedPrec);

        cout << "  Trials:       " << trials << "\n"
             << "  Fail thresh:  " << failThreshold << "\n\n";
    }

    if (strategy == PlacementStrategy::GREEDY)
        initialiseGreedyCache(pattern, maxReadLength);

    // Matrix pre-filled with UINT16_MAX (pending sentinel) by dispatch section
    vector<uint16_t> matrix(rows * cols, UINT16_MAX);

    // ── Dispatch ALL cells to thread pool ────────────────────────────────
    // Workers pull cells from atomic counter — fully parallel, lock-free.
    // Each worker peeks left in its row for warm start bounds:
    //   - UINT16_MAX = pending (skip)
    //   - v <= maxS  = viable (use as hi upper bound)
    //   - v > maxS   = unviable sentinel maxS+1 (use maxS as lo lower bound)
    //
    // No shared mutable state — workers only write their own cell.
    // uint16_t reads/writes are atomic on x86.

    // Pre-fill with UINT16_MAX (pending sentinel)
    fill(matrix.begin(), matrix.end(), UINT16_MAX);

    struct WorkItem {
        int r, c;
    };
    vector<WorkItem> work;
    work.reserve(rows * cols);
    for (int r = 0; r < rows; r++)
        for (int c = 0; c < cols; c++)
            work.push_back({r, c});

    atomic<int> nextCell(0);
    atomic<int> completedCells(0);
    mutex printMutex;

    const int numThreads = max(1u, thread::hardware_concurrency());
    const int minViable = max(3, (GATE_COVERAGE + pattern.weight - 1) / pattern.weight);

    vector<thread> workers;
    for (int t = 0; t < numThreads; t++) {
        workers.emplace_back([&]() {
            while (true) {
                int idx = nextCell++;
                if (idx >= (int)work.size())
                    break;

                const int r = work[idx].r;
                const int c = work[idx].c;
                const double errorRate = errorRates[r];
                const int maxSeeds = readLengths[c] - pattern.window + 1;

                // Peek left for warm start bounds
                int hiWarm = maxSeeds;  // default: full range
                int loWarm = minViable; // default: minimum

                for (int cp = c - 1; cp >= 0; cp--) {
                    uint16_t v = matrix[r * cols + cp];
                    if (v == UINT16_MAX)
                        continue; // pending — skip

                    int maxS = readLengths[cp] - pattern.window + 1;
                    if ((int)v <= maxS) {
                        // Viable cell — use budget as upper bound
                        hiWarm = min(hiWarm, (int)v);
                    }
                    else {
                        // Unviable cell — use maxS as lower bound
                        loWarm = max(loWarm, maxS);
                    }
                    break; // stop at first non-pending cell
                }

                int budget = MonteCarlo::findMinimumSeeds(pattern, readLengths[c], errorRate,
                                                          targetProbability, strategy, errorType,
                                                          trials, failThreshold, hiWarm);

                // Store result — replaces UINT16_MAX pending sentinel
                // budget <= maxSeeds: achievable
                // budget == maxSeeds+1: unachievable sentinel
                matrix[r * cols + c] = (uint16_t)budget;

                int done = ++completedCells;
                if ((done % 1000) == 0 || done == (int)work.size()) {
                    lock_guard<mutex> lk(printMutex);
                    cerr << "\r[" << done << "/" << work.size() << "] " << fixed << setprecision(1)
                         << (100.0 * done / work.size()) << "%" << flush;
                }
            }
        });
    }
    for (auto& th : workers)
        th.join();
    cerr << "\n";

    // Note: matrix contains maxSeeds+1 sentinels for unachievable cells.
    // These are clamped to maxSeeds when writing to bin for skippy compatibility.

    // ── Statistics ────────────────────────────────────────────────────────
    double sum = 0, frob = 0;
    uint16_t maxCell = 0;
    int achievable = 0;
    for (int r = 0; r < rows; r++)
        for (int c = 0; c < cols; c++) {
            int maxS = readLengths[c] - pattern.window + 1;
            uint16_t v = matrix[r * cols + c];
            if ((int)v <= maxS) {
                sum += v;
                frob += (double)v * v;
                if (v > maxCell)
                    maxCell = v;
                achievable++;
            }
        }
    double mean = achievable > 0 ? sum / achievable : 0;
    frob = sqrt(frob);

    // ── Build and save SeedMatrix ─────────────────────────────────────────
    SeedMatrix sm;
    // Clamp sentinels to maxSeeds before storing in SeedMatrix
    // so skippy always gets feasible seed counts
    vector<uint16_t> clampedMatrix(matrix);
    for (int r = 0; r < rows; r++)
        for (int c = 0; c < cols; c++) {
            int idx = r * cols + c;
            int maxS = readLengths[c] - pattern.window + 1;
            if ((int)clampedMatrix[idx] > maxS)
                clampedMatrix[idx] = (uint16_t)maxS;
        }

    sm.setReadLengths(readLengths);
    sm.setErrorRates(errorRates);
    sm.setMatrix(clampedMatrix);

    auto& hdr = sm.getHeader();
    strncpy(hdr.pattern, pattern.pattern_str.c_str(), sizeof(hdr.pattern) - 1);
    hdr.pattern[sizeof(hdr.pattern) - 1] = '\0';
    hdr.version = 3;
    hdr.targetProbability = targetProbability;
    hdr.targetCoverage = GATE_COVERAGE;
    hdr.seedSpan = pattern.window;
    hdr.informativeBases = pattern.weight;

    string binName = baseName + suffix + ".bin";
    if (!sm.save(binName)) {
        cerr << "ERROR: Failed to write " << binName << "\n";
        return false;
    }

    string csvName = baseName + suffix + ".csv";
    ofstream csv(csvName);
    if (csv.is_open()) {
        csv << "ErrorRate";
        for (int w : readLengths)
            csv << "," << w;
        csv << "\n";
        for (int r = 0; r < rows; r++) {
            csv << errorRates[r];
            for (int c = 0; c < cols; c++)
                csv << "," << matrix[r * cols + c];
            csv << "\n";
        }
    }

    {
        ios_base::fmtflags savedFlags = cout.flags();
        streamsize savedPrec = cout.precision();
        cout << "Pattern:     " << pattern.pattern_str << "\n"
             << "Error type:  " << (errorType == ErrorType::INDEL ? "indels" : "subs") << "\n"
             << "Achievable:  " << achievable << "/" << (rows * cols) << " cells\n"
             << fixed << setprecision(1) << "Mean budget: " << mean << "\n"
             << "Frobenius:   " << frob << "\n"
             << "Max budget:  " << maxCell << "\n"
             << "Saved:       " << binName << "\n\n";
        cout.flags(savedFlags);
        cout.precision(savedPrec);
    }

    // ── Build truncation matrix by reading back the CSV ─────────────────
    // (mirrors the standalone test_truncation logic which is verified correct)
    {
        // Re-read the CSV we just wrote
        ifstream csvIn(csvName);
        vector<int> csvReadLengths;
        vector<double> csvErrorRates;
        vector<uint16_t> csvMatrix;

        string line;
        bool firstLine = true;
        while (getline(csvIn, line)) {
            stringstream ss(line);
            string token;
            vector<string> tokens;
            while (getline(ss, token, ','))
                tokens.push_back(token);
            if (firstLine) {
                for (size_t i = 1; i < tokens.size(); i++)
                    csvReadLengths.push_back(stoi(tokens[i]));
                firstLine = false;
                continue;
            }
            csvErrorRates.push_back(stod(tokens[0]));
            for (size_t i = 1; i < tokens.size(); i++)
                csvMatrix.push_back((uint16_t)stoi(tokens[i]));
        }

        const int csvRows = (int)csvErrorRates.size();
        const int csvCols = (int)csvReadLengths.size();
        const int csvSpan = pattern.window;

        vector<uint16_t> truncLen(csvRows * csvCols, 0);
        vector<uint16_t> truncSeeds(csvRows * csvCols, 0);
        vector<uint16_t> truncErrorRate(csvRows * csvCols, 0);

        for (int r = 0; r < csvRows; r++) {
            int dominantLen = -1;
            int dominantSeeds = INT_MAX;

            for (int c = 0; c < csvCols; c++) {
                int idx = r * csvCols + c;
                int maxS = csvReadLengths[c] - csvSpan + 1;
                int v = (int)csvMatrix[idx];

                if (v > maxS) // unachievable sentinel
                {
                    truncLen[idx] = (uint16_t)csvReadLengths[c];
                    truncSeeds[idx] = (uint16_t)maxS;
                    truncErrorRate[idx] = (uint16_t)csvRows;
                    continue;
                }

                if (v < dominantSeeds) {
                    dominantSeeds = v;
                    dominantLen = csvReadLengths[c];
                }

                truncLen[idx] = (uint16_t)dominantLen;
                truncSeeds[idx] = (uint16_t)dominantSeeds;
                truncErrorRate[idx] = (uint16_t)r;
            }
        }

        // Store in SeedMatrix and save updated bin
        sm.setTruncLengths(truncLen);
        sm.setTruncSeeds(truncSeeds);
        sm.setTruncErrorRates(truncErrorRate);

        if (!sm.save(binName))
            cerr << "ERROR: Failed to write truncation to " << binName << "\n";
        else
            cout << "Updated " << binName << " with truncation matrix\n";

        // Save truncation CSVs
        auto saveTruncCSV = [&](const string& path, const vector<uint16_t>& data) {
            ofstream f(path);
            f << "ErrorRate";
            for (int w : csvReadLengths)
                f << "," << w;
            f << "\n";
            for (int r = 0; r < csvRows; r++) {
                f << csvErrorRates[r];
                for (int c = 0; c < csvCols; c++)
                    f << "," << data[r * csvCols + c];
                f << "\n";
            }
            cout << "Saved " << path << "\n";
        };

        saveTruncCSV(baseName + suffix + "_trunc_len.csv", truncLen);
        saveTruncCSV(baseName + suffix + "_trunc_seeds.csv", truncSeeds);
        saveTruncCSV(baseName + suffix + "_trunc_errrate.csv", truncErrorRate);

        // Simple contract report
        int satCells = 0, truncBenefitCells = 0;
        double sumSeeds = 0, sumSq = 0;
        for (int r = 0; r < csvRows; r++)
            for (int c = 0; c < csvCols; c++) {
                int idx = r * csvCols + c;
                if ((int)truncErrorRate[idx] == r) {
                    satCells++;
                    int n = (int)truncSeeds[idx];
                    sumSeeds += n;
                    sumSq += (double)n * n;
                    if ((int)truncLen[idx] < csvReadLengths[c])
                        truncBenefitCells++;
                }
            }

        double mean = satCells > 0 ? sumSeeds / satCells : 0;
        double rms = satCells > 0 ? sqrt(sumSq / satCells) : 0;

        cout << "\n====================================\n"
             << "Contract: " << pattern.pattern_str << " " << suffix.substr(1) << "\n"
             << "  Grid:         " << csvRows << " x " << csvCols << "\n"
             << "  Satisfiable:  " << satCells << "/" << csvRows * csvCols << "\n"
             << "  Mean seeds:   " << mean << "\n"
             << "  RMS seeds:    " << rms << "\n"
             << "  Trunc saving: " << truncBenefitCells << " cells\n"
             << "====================================\n\n";
    }

    return true;
}

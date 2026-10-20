// dump_trunc_matrix.cpp
// Loads a version 2 .bin file and dumps the truncation matrices to CSV.
//
// Usage: ./dump_trunc_matrix <binfile> [trunc_len.csv] [trunc_seeds.csv]
//
// Example:
//   ./dump_trunc_matrix GRCh38_greedy.bin trunc_len.csv trunc_seeds.csv

#include <fstream>
#include <iostream>
#include <string>

#include "SeedMatrix.hpp"

using namespace std;

static void saveCSV(const string& path, const vector<uint16_t>& data,
                    const vector<int>& readLengths, const vector<double>& errorRates)
{
    const int rows = (int)errorRates.size();
    const int cols = (int)readLengths.size();

    ofstream f(path);
    if (!f) {
        cerr << "Cannot open " << path << "\n";
        return;
    }

    f << "ErrorRate";
    for (int w : readLengths)
        f << "," << w;
    f << "\n";

    for (int r = 0; r < rows; r++) {
        f << errorRates[r];
        for (int c = 0; c < cols; c++)
            f << "," << data[r * cols + c];
        f << "\n";
    }
    cout << "Saved " << path << "\n";
}

int main(int argc, char* argv[])
{
    if (argc < 2) {
        cerr << "Usage: " << argv[0] << " <binfile> [trunc_len.csv] [trunc_seeds.csv]\n";
        return 1;
    }

    string binFile = argv[1];
    string lenCsvFile = (argc > 2) ? argv[2] : binFile + "_trunc_len.csv";
    string seedsCsvFile = (argc > 3) ? argv[3] : binFile + "_trunc_seeds.csv";

    SeedMatrix sm;
    if (!sm.load(binFile)) {
        cerr << "Failed to load " << binFile << "\n";
        return 1;
    }

    sm.printSummary();

    if (!sm.hasTruncMatrix()) {
        cerr << "No truncation matrix in " << binFile << " (version " << sm.getHeader().version
             << ")\n";
        return 1;
    }

    // Access truncation data via queryTrunc across the full grid
    const auto& readLengths = sm.getReadLengths();
    const auto& errorRates = sm.getErrorRates();
    const int rows = (int)errorRates.size();
    const int cols = (int)readLengths.size();

    vector<uint16_t> truncLen(rows * cols);
    vector<uint16_t> truncSeeds(rows * cols);

    for (int r = 0; r < rows; r++) {
        for (int c = 0; c < cols; c++) {
            auto ts = sm.queryTrunc(errorRates[r], readLengths[c]);
            truncLen[r * cols + c] = (uint16_t)ts.truncLength;
            truncSeeds[r * cols + c] = (uint16_t)ts.seedCount;
        }
    }

    saveCSV(lenCsvFile, truncLen, readLengths, errorRates);
    saveCSV(seedsCsvFile, truncSeeds, readLengths, errorRates);

    return 0;
}

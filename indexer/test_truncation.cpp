// test_truncation.cpp
// Reads GRCh38_greedy.csv, builds truncation matrix, writes trunc_seeds CSV.
// Standalone test of buildTruncationMatrix only.

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <climits>
#include <cmath>
#include <algorithm>
#include <cstdint>

using namespace std;

int main(int argc, char* argv[])
{
    if (argc < 4)
    {
        cerr << "Usage: test_truncation <input.csv> <window> <output.csv>\n";
        cerr << "  e.g. test_truncation GRCh38_greedy.csv 24 trunc_seeds_standalone.csv\n";
        return 1;
    }

    const string inFile  = argv[1];
    const int    span    = atoi(argv[2]);
    const string outFile = argv[3];

    // ── Read CSV ──────────────────────────────────────────────────────────
    ifstream f(inFile);
    if (!f) { cerr << "Cannot open " << inFile << "\n"; return 1; }

    vector<double>   errorRates;
    vector<int>      readLengths;
    vector<uint16_t> matrix;

    string line;
    bool firstLine = true;
    while (getline(f, line))
    {
        stringstream ss(line);
        string token;
        vector<string> tokens;
        while (getline(ss, token, ','))
            tokens.push_back(token);

        if (firstLine)
        {
            // Header row: ErrorRate, 60, 61, 62, ...
            for (size_t i = 1; i < tokens.size(); i++)
                readLengths.push_back(stoi(tokens[i]));
            firstLine = false;
            continue;
        }

        errorRates.push_back(stod(tokens[0]));
        for (size_t i = 1; i < tokens.size(); i++)
            matrix.push_back((uint16_t)stoi(tokens[i]));
    }

    const int rows = (int)errorRates.size();
    const int cols = (int)readLengths.size();

    cout << "Read " << rows << " error rates x " << cols << " read lengths\n";

    // ── Build truncation matrix ───────────────────────────────────────────
    vector<uint16_t> truncLen(rows * cols, 0);
    vector<uint16_t> truncSeeds(rows * cols, 0);

    for (int r = 0; r < rows; r++)
    {
        int dominantLen   = -1;
        int dominantSeeds = INT_MAX;

        for (int c = 0; c < cols; c++)
        {
            int idx  = r * cols + c;
            int maxS = readLengths[c] - span + 1;

            if ((int)matrix[idx] > maxS)  // unachievable sentinel
            {
                truncLen  [idx] = (uint16_t)readLengths[c];
                truncSeeds[idx] = (uint16_t)maxS;
                continue;
            }

            // Achievable
            if ((int)matrix[idx] < dominantSeeds)
            {
                dominantSeeds = (int)matrix[idx];
                dominantLen   = readLengths[c];
            }

            truncLen  [idx] = (uint16_t)dominantLen;
            truncSeeds[idx] = (uint16_t)dominantSeeds;
        }
    }

    // ── Write CSV ─────────────────────────────────────────────────────────
    ofstream out(outFile);
    out << "ErrorRate";
    for (int w : readLengths) out << "," << w;
    out << "\n";
    for (int r = 0; r < rows; r++)
    {
        out << errorRates[r];
        for (int c = 0; c < cols; c++)
            out << "," << truncSeeds[r * cols + c];
        out << "\n";
    }

    cout << "Written to " << outFile << "\n";
    return 0;
}

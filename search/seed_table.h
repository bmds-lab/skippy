// SeedTable: computes minimum seeds needed to pass the coverage
// gate with target probability, for a grid of (error_rate, read_length).
// Computed once at startup using Monte Carlo simulation.

#pragma once
#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iostream>
#include <mutex>
#include <random>
#include <thread>
#include <vector>

class SeedTable {
public:
    static constexpr int WINDOW = 21;
    static constexpr int WEIGHT = 16;
    static constexpr int USED[16] = {0, 1, 3, 4, 5, 7, 8, 10, 11, 12, 13, 15, 16, 18, 19, 20};
    static constexpr int N_SIM = 300; // per (K, L, e) point — fast but sufficient

    static constexpr int READ_LENGTHS[13] = {60,  80,  100, 120, 150, 200, 250,
                                             300, 400, 500, 600, 800, 1000};
    static constexpr float ERROR_RATES[17] = {1,  2,  3,  4,  5,  6,  7,  8, 9,
                                              10, 12, 14, 15, 17, 20, 22, 25};
    static constexpr int N_LENGTHS = 13;
    static constexpr int N_ERRORS = 17;

    // Compute the table. Call once at startup.
    // gate:        skippy gate parameter (e.g. 28)
    // targetProb:  probability threshold (e.g. 0.99)
    // fudge:       safety multiplier (e.g. 1.2 to be conservative)
    // nThreads:    parallelism (use skippy's THREADS global)

    void clear()
    {
        table_.assign(N_ERRORS, std::vector<int>(N_LENGTHS, -1));
    }

    void setCell(int errorIndex, int lengthIndex, int value)
    {
        table_[errorIndex][lengthIndex] = value;
    }

    bool save(const std::string& filename) const
    {
        std::ofstream out(filename, std::ios::binary);

        if (!out)
            return false;

        int rows = table_.size();
        int cols = rows ? table_[0].size() : 0;

        out.write((char*)&rows, sizeof(rows));
        out.write((char*)&cols, sizeof(cols));

        for (int r = 0; r < rows; ++r)
            out.write((char*)table_[r].data(), cols * sizeof(int));

        return true;
    }

    bool load(const std::string& filename)
    {
        std::ifstream in(filename, std::ios::binary);

        if (!in)
            return false;

        int rows, cols;

        in.read((char*)&rows, sizeof(rows));
        in.read((char*)&cols, sizeof(cols));

        if (!in)
            return false;

        table_.assign(rows, std::vector<int>(cols));

        for (int r = 0; r < rows; ++r)
            in.read((char*)table_[r].data(), cols * sizeof(int));

        return (bool)in;
    }
    void compute(uint32_t gate, float targetProb, float fudge = 1.0f, int nThreads = 1)
    {
        gate_ = gate;
        targetProb_ = targetProb;
        fudge_ = fudge;

        auto t0 = std::chrono::steady_clock::now();

        table_.assign(N_ERRORS, std::vector<int>(N_LENGTHS, -1));

        // Launch one thread per error-rate row
        int nT = std::min(nThreads, N_ERRORS);
        std::vector<std::thread> threads;
        std::mutex mtx;
        int nextRow = 0;

        for (int t = 0; t < nT; t++) {
            threads.emplace_back([&, t]() {
                std::mt19937 rng(42 + t);
                while (true) {
                    int ei;
                    {
                        std::lock_guard<std::mutex> lk(mtx);
                        ei = nextRow++;
                    }
                    if (ei >= N_ERRORS)
                        break;

                    float errRate = ERROR_RATES[ei] / 100.0f;
                    for (int li = 0; li < N_LENGTHS; li++) {
                        int L = READ_LENGTHS[li];
                        table_[ei][li] = findMinSeeds(L, errRate, gate, targetProb, fudge, rng);
                    }
                }
            });
        }
        for (auto& th : threads)
            th.join();

        auto t1 = std::chrono::steady_clock::now();
        double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
        std::cout << "SeedTable computed in " << ms << " ms"
                  << " (" << nT << " threads)\n";
    }

    // Given error rate (percent) and read length, return seeds needed.
    uint32_t lookup(float errorPct, uint32_t readLen) const
    {
        if (table_.empty())
            return 80;

        // Clamp to grid
        float e = std::max(ERROR_RATES[0], std::min(ERROR_RATES[N_ERRORS - 1], errorPct));
        float L = std::max((float)READ_LENGTHS[0],
                           std::min((float)READ_LENGTHS[N_LENGTHS - 1], (float)readLen));

        // Bracket error rate
        int ei = 0;
        while (ei < N_ERRORS - 2 && ERROR_RATES[ei + 1] < e)
            ei++;
        float te = (e - ERROR_RATES[ei]) / (ERROR_RATES[ei + 1] - ERROR_RATES[ei]);

        // Bracket read length
        int li = 0;
        while (li < N_LENGTHS - 2 && READ_LENGTHS[li + 1] < (int)L)
            li++;
        float tl = (L - READ_LENGTHS[li]) / (float)(READ_LENGTHS[li + 1] - READ_LENGTHS[li]);

        // Bilinear interpolation — use maxPossible for unachievable cells
        float v00 = cellVal(ei, li);
        float v10 = cellVal(ei + 1, li);
        float v01 = cellVal(ei, li + 1);
        float v11 = cellVal(ei + 1, li + 1);

        float v =
            (1 - te) * (1 - tl) * v00 + te * (1 - tl) * v10 + (1 - te) * tl * v01 + te * tl * v11;

        uint32_t maxPoss = (readLen >= (uint32_t)WINDOW) ? (readLen - WINDOW + 1) : 0;
        return std::min((uint32_t)std::ceil(v), maxPoss);
    }

    void print() const
    {
        std::cout << "SeedTable (gate=" << gate_ << " P=" << targetProb_ << " fudge=" << fudge_
                  << ")\n";
        std::cout << "err%";
        for (int li = 0; li < N_LENGTHS; li++)
            std::cout << "\t" << READ_LENGTHS[li];
        std::cout << "\n";
        for (int ei = 0; ei < N_ERRORS; ei++) {
            std::cout << ERROR_RATES[ei] << "%";
            for (int li = 0; li < N_LENGTHS; li++) {
                int v = table_[ei][li];
                if (v < 0)
                    std::cout << "\t—";
                else
                    std::cout << "\t" << v;
            }
            std::cout << "\n";
        }
    }

private:
    uint32_t gate_ = 28;
    float targetProb_ = 0.99f;
    float fudge_ = 1.0f;
    std::vector<std::vector<int>> table_;

    float cellVal(int ei, int li) const
    {
        int v = table_[ei][li];
        return (v < 0) ? (float)(READ_LENGTHS[li] - WINDOW + 1) : (float)v;
    }

    float simProb(int L, int K, float errRate, uint32_t gate, std::mt19937& rng) const
    {
        int maxPoss = L - WINDOW + 1;
        if (maxPoss <= 0 || K <= 0)
            return 0.0f;
        K = std::min(K, maxPoss);

        std::vector<int> seedPos(K);
        for (int i = 0; i < K; i++)
            seedPos[i] = (K == 1) ? 0 : (int)std::round(i * (maxPoss - 1) / (double)(K - 1));

        std::uniform_real_distribution<float> uni(0.0f, 1.0f);
        std::vector<uint8_t> errors(L), bitmap(L);
        int pass = 0;

        for (int s = 0; s < N_SIM; s++) {
            for (int b = 0; b < L; b++)
                errors[b] = uni(rng) < errRate ? 1 : 0;

            std::fill(bitmap.begin(), bitmap.end(), 0);
            uint32_t coverage = 0;

            for (int si = 0; si < K; si++) {
                int seed = seedPos[si];
                bool intact = true;
                for (int k = 0; k < WEIGHT; k++)
                    if (errors[seed + USED[k]]) {
                        intact = false;
                        break;
                    }
                if (!intact)
                    continue;
                for (int k = 0; k < WEIGHT; k++) {
                    int bit = seed + USED[k];
                    if (!bitmap[bit]) {
                        bitmap[bit] = 1;
                        coverage++;
                    }
                }
            }
            if (coverage >= gate)
                pass++;
        }
        return (float)pass / N_SIM;
    }

    int findMinSeeds(int L, float errRate, uint32_t gate, float targetProb, float fudge,
                     std::mt19937& rng) const
    {
        int maxPoss = L - WINDOW + 1;
        if (maxPoss <= 0)
            return -1;
        if (simProb(L, maxPoss, errRate, gate, rng) < targetProb)
            return -1;

        int step = std::max(1, maxPoss / 30);
        int prev = 0;
        for (int n = step; n <= maxPoss; n += step) {
            if (simProb(L, n, errRate, gate, rng) >= targetProb) {
                for (int nn = prev + 1; nn < n; nn++)
                    if (simProb(L, nn, errRate, gate, rng) >= targetProb)
                        return (int)std::ceil(nn * fudge);
                return (int)std::ceil(n * fudge);
            }
            prev = n;
        }
        return (int)std::ceil(maxPoss * fudge);
    }
};

constexpr int SeedTable::USED[];
constexpr int SeedTable::READ_LENGTHS[];
constexpr float SeedTable::ERROR_RATES[];

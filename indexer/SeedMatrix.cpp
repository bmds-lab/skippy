// SeedMatrix.cpp — shared implementation (indexer and search)
#include "SeedMatrix.hpp"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <fstream>
#include <iostream>
using namespace std;

SeedMatrix::SeedMatrix()
{
    clear();
}

void SeedMatrix::clear()
{
    memset(&header_, 0, sizeof(header_));
    header_.version = 3;
    readLengths_.clear();
    errorRates_.clear();
    matrix_.clear();
    truncLengths_.clear();
    truncSeeds_.clear();
    truncErrorRates_.clear();
}

const vector<int>& SeedMatrix::getReadLengths() const
{
    return readLengths_;
}
const vector<double>& SeedMatrix::getErrorRates() const
{
    return errorRates_;
}
const vector<uint16_t>& SeedMatrix::getMatrix() const
{
    return matrix_;
}
SeedMatrix::Header& SeedMatrix::getHeader()
{
    return header_;
}
const SeedMatrix::Header& SeedMatrix::getHeader() const
{
    return header_;
}

void SeedMatrix::setReadLengths(const vector<int>& v)
{
    readLengths_ = v;
    header_.numReadLengths = (uint32_t)v.size();
}
void SeedMatrix::setErrorRates(const vector<double>& v)
{
    errorRates_ = v;
    header_.numErrorRates = (uint32_t)v.size();
}
void SeedMatrix::setMatrix(const vector<uint16_t>& v)
{
    matrix_ = v;
}
void SeedMatrix::setTruncLengths(const vector<uint16_t>& v)
{
    truncLengths_ = v;
}
void SeedMatrix::setTruncSeeds(const vector<uint16_t>& v)
{
    truncSeeds_ = v;
}
void SeedMatrix::setTruncErrorRates(const vector<uint16_t>& v)
{
    truncErrorRates_ = v;
}

bool SeedMatrix::hasTruncMatrix() const
{
    return header_.version >= 3 && !truncLengths_.empty() && !truncSeeds_.empty() &&
           !truncErrorRates_.empty();
}

bool SeedMatrix::save(const std::string& filename) const
{
    ofstream out(filename, ios::binary);
    if (!out)
        return false;
    Header hdr = header_;
    hdr.version = hasTruncMatrix() ? 3 : 1;
    out.write((const char*)&hdr, sizeof(hdr));
    out.write((const char*)readLengths_.data(), readLengths_.size() * sizeof(int));
    out.write((const char*)errorRates_.data(), errorRates_.size() * sizeof(double));
    out.write((const char*)matrix_.data(), matrix_.size() * sizeof(uint16_t));
    if (hdr.version >= 3) {
        size_t n = header_.numReadLengths * header_.numErrorRates;
        out.write((const char*)truncLengths_.data(), n * sizeof(uint16_t));
        out.write((const char*)truncSeeds_.data(), n * sizeof(uint16_t));
        out.write((const char*)truncErrorRates_.data(), n * sizeof(uint16_t));
    }
    out.close();
    return static_cast<bool>(out);
}

bool SeedMatrix::load(const string& filename)
{
    clear();
    ifstream in(filename, ios::binary);
    if (!in)
        return false;
    in.read((char*)&header_, sizeof(header_));
    readLengths_.resize(header_.numReadLengths);
    errorRates_.resize(header_.numErrorRates);
    matrix_.resize(header_.numReadLengths * header_.numErrorRates);
    in.read((char*)readLengths_.data(), readLengths_.size() * sizeof(int));
    in.read((char*)errorRates_.data(), errorRates_.size() * sizeof(double));
    in.read((char*)matrix_.data(), matrix_.size() * sizeof(uint16_t));
    if (header_.version >= 3) {
        size_t n = header_.numReadLengths * header_.numErrorRates;
        truncLengths_.resize(n);
        truncSeeds_.resize(n);
        truncErrorRates_.resize(n);
        in.read((char*)truncLengths_.data(), n * sizeof(uint16_t));
        in.read((char*)truncSeeds_.data(), n * sizeof(uint16_t));
        in.read((char*)truncErrorRates_.data(), n * sizeof(uint16_t));
    }
    return true;
}

int SeedMatrix::findLowerReadIndex(int w) const
{
    if (w <= readLengths_.front())
        return 0;
    if (w >= readLengths_.back())
        return (int)readLengths_.size() - 2;
    for (size_t i = 0; i + 1 < readLengths_.size(); i++)
        if (readLengths_[i] <= w && w <= readLengths_[i + 1])
            return (int)i;
    return (int)readLengths_.size() - 2;
}

int SeedMatrix::findLowerErrorIndex(double e) const
{
    if (e <= errorRates_.front())
        return 0;
    if (e >= errorRates_.back())
        return (int)errorRates_.size() - 2;
    for (size_t i = 0; i + 1 < errorRates_.size(); i++)
        if (errorRates_[i] <= e && e <= errorRates_[i + 1])
            return (int)i;
    return (int)errorRates_.size() - 2;
}

int SeedMatrix::findNearestReadIndex(int w) const
{
    int best = 0, bestDist = abs(w - readLengths_[0]);
    for (size_t i = 1; i < readLengths_.size(); i++) {
        int d = abs(w - readLengths_[i]);
        if (d < bestDist) {
            bestDist = d;
            best = (int)i;
        }
    }
    return best;
}

int SeedMatrix::findNearestErrorIndex(double e) const
{
    int best = 0;
    double bestDist = fabs(e - errorRates_[0]);
    for (size_t i = 1; i < errorRates_.size(); i++) {
        double d = fabs(e - errorRates_[i]);
        if (d < bestDist) {
            bestDist = d;
            best = (int)i;
        }
    }
    return best;
}

double SeedMatrix::bilinearInterpolate(double e, int w) const
{
    int r = findLowerErrorIndex(e);
    int c = findLowerReadIndex(w);
    double e0 = errorRates_[r], e1 = errorRates_[r + 1];
    double w0 = readLengths_[c], w1 = readLengths_[c + 1];
    auto value = [&](int rr, int cc) { return (double)matrix_[rr * header_.numReadLengths + cc]; };
    double tx = (w - w0) / (w1 - w0);
    double ty = (e - e0) / (e1 - e0);
    double top = value(r, c) * (1.0 - tx) + value(r, c + 1) * tx;
    double bot = value(r + 1, c) * (1.0 - tx) + value(r + 1, c + 1) * tx;
    return top * (1.0 - ty) + bot * ty;
}

int SeedMatrix::query(double errorRate, int readLength) const
{
    return (int)ceil(bilinearInterpolate(errorRate, readLength));
}

TruncatedSeed SeedMatrix::queryTrunc(double errorRate, int readLength) const
{
    if (!hasTruncMatrix())
        return {readLength, query(errorRate, readLength)};

    int r = findNearestErrorIndex(errorRate);
    int c = findNearestReadIndex(readLength);
    int idx = r * (int)header_.numReadLengths + c;

    int tlen = (int)truncLengths_[idx];
    int tseeds = (int)truncSeeds_[idx];

    if (tlen > readLength)
        tlen = readLength;

    return {tlen, tseeds};
}

void SeedMatrix::printSummary() const
{
    cout << "Pattern: " << header_.pattern << "\n"
         << "Version: " << header_.version << "\n"
         << "Truncation: " << (hasTruncMatrix() ? "yes" : "no") << "\n";
}

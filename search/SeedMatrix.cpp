#include "SeedMatrix.hpp"

#include <fstream>
#include <iostream>
#include <cstring>
#include <cmath>
#include <algorithm>

using namespace std;


/*=========================================================
    Constructor
=========================================================*/

SeedMatrix::SeedMatrix()
{
    clear();
}


/*=========================================================
    Clear
=========================================================*/

void SeedMatrix::clear()
{
    memset(&header_, 0, sizeof(header_));

    header_.version = 1;

    readLengths_.clear();
    errorRates_.clear();
    matrix_.clear();
}


/*=========================================================
    Save
=========================================================*/

bool SeedMatrix::save(
    const string& filename) const
{
    ofstream out(
        filename,
        ios::binary);

    if (!out)
    {
        return false;
    }

    out.write(
        (const char*)&header_,
        sizeof(header_));

    out.write(
        (const char*)readLengths_.data(),
        readLengths_.size()
            * sizeof(int));

    out.write(
        (const char*)errorRates_.data(),
        errorRates_.size()
            * sizeof(double));

    out.write(
        (const char*)matrix_.data(),
        matrix_.size()
            * sizeof(uint16_t));

    return true;
}


/*=========================================================
    Load
=========================================================*/

bool SeedMatrix::load(
    const string& filename)
{
    clear();

    ifstream in(
        filename,
        ios::binary);

    if (!in)
    {
        return false;
    }

    in.read(
        (char*)&header_,
        sizeof(header_));
/*
cerr << "pattern          = " << header_.pattern << endl;
cerr << "version          = " << header_.version << endl;
cerr << "targetProb       = " << header_.targetProbability << endl;
cerr << "targetCoverage   = " << header_.targetCoverage << endl;
cerr << "seedSpan         = " << header_.seedSpan << endl;
cerr << "informativeBases = " << header_.informativeBases << endl;
cerr << "numReadLengths   = " << header_.numReadLengths << endl;
cerr << "numErrorRates    = " << header_.numErrorRates << endl;
*/

    readLengths_.resize(
        header_.numReadLengths);

    errorRates_.resize(
        header_.numErrorRates);

    matrix_.resize(
        header_.numReadLengths
        *
        header_.numErrorRates);

    in.read(
        (char*)readLengths_.data(),
        readLengths_.size()
            * sizeof(int));

    in.read(
        (char*)errorRates_.data(),
        errorRates_.size()
            * sizeof(double));

    in.read(
        (char*)matrix_.data(),
        matrix_.size()
            * sizeof(uint16_t));
/*
cerr << "readLengths_.size() = "
     << readLengths_.size()
     << endl;

cerr << "errorRates_.size() = "
     << errorRates_.size()
     << endl;

cerr << "matrix_.size() = "
     << matrix_.size()
     << endl;

cerr << "first length = "
     << readLengths_.front()
     << " last length = "
     << readLengths_.back()
     << endl;

cerr << "first error = "
     << errorRates_.front()
     << " last error = "
     << errorRates_.back()
     << endl;
*/

    if (header_.version >= 2)
    {
        size_t n = header_.numReadLengths * header_.numErrorRates;
        truncLengths_.resize(n);
        truncSeeds_.resize(n);
        in.read((char*)truncLengths_.data(), n * sizeof(uint16_t));
        in.read((char*)truncSeeds_.data(),   n * sizeof(uint16_t));
    }

    return true;
}


/*=========================================================
    Accessors
=========================================================*/

const vector<int>&
SeedMatrix::getReadLengths() const
{
    return readLengths_;
}

const vector<double>&
SeedMatrix::getErrorRates() const
{
    return errorRates_;
}

const vector<uint16_t>&
SeedMatrix::getMatrix() const
{
    return matrix_;
}

SeedMatrix::Header&
SeedMatrix::getHeader()
{
    return header_;
}

const SeedMatrix::Header&
SeedMatrix::getHeader() const
{
    return header_;
}

void SeedMatrix::setReadLengths(
    const vector<int>& v)
{
    readLengths_ = v;

    header_.numReadLengths =
        (uint32_t)v.size();
}

void SeedMatrix::setErrorRates(
    const vector<double>& v)
{
    errorRates_ = v;

    header_.numErrorRates =
        (uint32_t)v.size();
}

void SeedMatrix::setMatrix(
    const vector<uint16_t>& v)
{
    matrix_ = v;
}


/*=========================================================
    Summary
=========================================================*/

void SeedMatrix::printSummary() const
{
    cout << "\nSeed Matrix\n";

    cout << "Pattern: "
         << header_.pattern
         << "\n";

    cout << "Probability: "
         << header_.targetProbability
         << "\n";

    cout << "Coverage: "
         << header_.targetCoverage
         << "\n";

    cout << "Read lengths: "
         << header_.numReadLengths
         << "\n";

    cout << "Error rates: "
         << header_.numErrorRates
         << "\n";
}


/*=========================================================
    Helper

    Largest index i such that

        readLengths_[i] <= w

=========================================================*/

int SeedMatrix::findLowerReadIndex(
    int w) const
{
    if (w <= readLengths_.front())
    {
        return 0;
    }

    if (w >= readLengths_.back())
    {
        return
            (int)readLengths_.size() - 2;
    }

    for (size_t i = 0;
         i + 1 < readLengths_.size();
         i++)
    {
        if (readLengths_[i] <= w
            &&
            w <= readLengths_[i+1])
        {
            return (int)i;
        }
    }

    return
        (int)readLengths_.size() - 2;
}


/*=========================================================
    Helper

    Largest index i such that

        errorRates_[i] <= e

=========================================================*/

int SeedMatrix::findLowerErrorIndex(
    double e) const
{
    if (e <= errorRates_.front())
    {
        return 0;
    }

    if (e >= errorRates_.back())
    {
        return
            (int)errorRates_.size() - 2;
    }

    for (size_t i = 0;
         i + 1 < errorRates_.size();
         i++)
    {
        if (errorRates_[i] <= e
            &&
            e <= errorRates_[i+1])
        {
            return (int)i;
        }
    }

    return
        (int)errorRates_.size() - 2;
}


/*=========================================================
    Bilinear interpolation
=========================================================*/

double SeedMatrix::bilinearInterpolate(
    double e,
    int w) const
{
    int r =
        findLowerErrorIndex(e);

    int c =
        findLowerReadIndex(w);

    double e0 =
        errorRates_[r];

    double e1 =
        errorRates_[r+1];

    double w0 =
        readLengths_[c];

    double w1 =
        readLengths_[c+1];

    auto value =
        [&](int rr,
            int cc)
        {
            return
            (double)
            matrix_[
                rr
                *
                header_.numReadLengths
                +
                cc
            ];
        };

    double q11 =
        value(r,c);

    double q12 =
        value(r,c+1);

    double q21 =
        value(r+1,c);

    double q22 =
        value(r+1,c+1);

    double tx =
        (w - w0)
        /
        (w1 - w0);

    double ty =
        (e - e0)
        /
        (e1 - e0);

    double top =
        q11 * (1.0 - tx)
        +
        q12 * tx;

    double bottom =
        q21 * (1.0 - tx)
        +
        q22 * tx;

    return
        top * (1.0 - ty)
        +
        bottom * ty;
}


/*=========================================================
    Query
=========================================================*/

int SeedMatrix::query(
    double errorRate,
    int readLength) const
{
    double v =
        bilinearInterpolate(
            errorRate,
            readLength);

    return
        (int)ceil(v);
}

// ── Truncation matrix support (version 2 .bin) ───────────────────────────────

void SeedMatrix::setTruncLengths(const std::vector<uint16_t>& v) { truncLengths_ = v; }
void SeedMatrix::setTruncSeeds(const std::vector<uint16_t>& v)   { truncSeeds_   = v; }

bool SeedMatrix::hasTruncMatrix() const
{
    return header_.version >= 2
           && !truncLengths_.empty()
           && !truncSeeds_.empty();
}

int SeedMatrix::findNearestReadIndex(int w) const
{
    int best = 0, bestDist = abs(w - readLengths_[0]);
    for (size_t i = 1; i < readLengths_.size(); i++)
    {
        int d = abs(w - readLengths_[i]);
        if (d < bestDist) { bestDist = d; best = (int)i; }
    }
    return best;
}

int SeedMatrix::findNearestErrorIndex(double e) const
{
    int    best     = 0;
    double bestDist = fabs(e - errorRates_[0]);
    for (size_t i = 1; i < errorRates_.size(); i++)
    {
        double d = fabs(e - errorRates_[i]);
        if (d < bestDist) { bestDist = d; best = (int)i; }
    }
    return best;
}

TruncatedSeed SeedMatrix::queryTrunc(double errorRate, int readLength) const
{
    if (!hasTruncMatrix())
        return { readLength, query(errorRate, readLength) };

    int r   = findNearestErrorIndex(errorRate);
    int c   = findNearestReadIndex(readLength);
    int idx = r * (int)header_.numReadLengths + c;

    int tlen   = (int)truncLengths_[idx];
    int tseeds = (int)truncSeeds_[idx];

    if (tlen > readLength) tlen = readLength;

    // Clamp seeds to what's actually available for the effective length
    int maxAvailable = (tlen >= header_.seedSpan)
        ? (tlen - (int)header_.seedSpan + 1)
        : 0;
    if (tseeds > maxAvailable) tseeds = maxAvailable;

    return { tlen, tseeds };
}

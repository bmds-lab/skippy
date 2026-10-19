// SeedMatrix.hpp — shared between indexer and search
#ifndef SEED_MATRIX_HPP
#define SEED_MATRIX_HPP

#include <string>
#include <vector>
#include <cstdint>

struct TruncatedSeed
{
    int truncLength;
    int seedCount;
};

class SeedMatrix
{
public:

    struct Header
    {
        char     pattern[64];
        uint32_t version;
        double   targetProbability;
        uint32_t targetCoverage;
        uint32_t seedSpan;
        uint32_t informativeBases;
        uint32_t numReadLengths;
        uint32_t numErrorRates;
    };

public:

    SeedMatrix();

    bool load(const std::string& filename);
    bool save(const std::string& filename) const;

    int           query(double errorRate, int readLength) const;
    TruncatedSeed queryTrunc(double errorRate, int readLength) const;
    bool          hasTruncMatrix() const;

    void clear();
    void printSummary() const;

    const std::vector<int>&      getReadLengths() const;
    const std::vector<double>&   getErrorRates()  const;
    const std::vector<uint16_t>& getMatrix()      const;

    Header&       getHeader();
    const Header& getHeader() const;

    void setReadLengths(const std::vector<int>& v);
    void setErrorRates(const std::vector<double>& v);
    void setMatrix(const std::vector<uint16_t>& v);
    void setTruncLengths(const std::vector<uint16_t>& v);
    void setTruncSeeds(const std::vector<uint16_t>& v);
    void setTruncErrorRates(const std::vector<uint16_t>& v);

private:

    int    findLowerReadIndex(int readLength)      const;
    int    findLowerErrorIndex(double errorRate)   const;
    int    findNearestReadIndex(int readLength)    const;
    int    findNearestErrorIndex(double errorRate) const;
    double bilinearInterpolate(double e, int w)    const;

private:

    Header                 header_;
    std::vector<int>       readLengths_;
    std::vector<double>    errorRates_;
    std::vector<uint16_t>  matrix_;
    std::vector<uint16_t>  truncLengths_;
    std::vector<uint16_t>  truncSeeds_;
    std::vector<uint16_t>  truncErrorRates_;  // error rate index certified for each cell
};

#endif

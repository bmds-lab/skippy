#include "SeedMatrix.hpp"
#include "generate_seed_matrix.hpp"
#include "indexGenome.hpp"

#include <cassert>
#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <string>
#include <vector>

static void checkPositions(size_t length, const std::string& first, const std::string& second = "")
{
    StripePattern pattern, other;
    if (!first.empty())
        assert(pattern.parse(first));
    if (!second.empty())
        assert(other.parse(second));
    std::string sequence;
    for (size_t i = 0; i < length; ++i)
        sequence += "ACGT"[i % 4];
    std::vector<protected_vector<uint32_t>> buckets(64);
    index_kmers(sequence.data(), length, buckets, 63, pattern, other);

    std::vector<unsigned> counts(length, 0);
    for (const auto& bucket : buckets) {
        for (uint32_t position : bucket) {
            assert(position < length);
            ++counts[position];
        }
    }
    const size_t window = first.empty() ? 32 : pattern.window;
    for (size_t i = 0; i < length; ++i) {
        unsigned expected = i + window <= length ? 1 : 0;
        if (!second.empty() && i + other.window <= length)
            ++expected;
        assert(counts[i] == expected);
    }
}

int main()
{
    for (size_t length : {0, 15, 16, 17, 31, 32, 64, 2200}) {
        checkPositions(length, "16");
        checkPositions(length, "8X8");
        checkPositions(length, "17");
        checkPositions(length, "");
        checkPositions(length, "16", "8XXXXX8");
        checkPositions(length, "8XXXXX8", "16");
    }

    char directory[] = "/tmp/skippy-index-test-XXXXXX";
    assert(mkdtemp(directory) != nullptr);
    const std::filesystem::path work(directory);
    SeedMatrix matrix;
    matrix.setReadLengths({60, 61});
    matrix.setErrorRates({0.0, 0.01});
    matrix.setMatrix({3, 4, 5, 6});
    assert(matrix.save((work / "matrix.bin").string()));
    SeedMatrix loaded;
    assert(loaded.load((work / "matrix.bin").string()));
    assert(loaded.getMatrix() == matrix.getMatrix());
    assert(!matrix.save("/dev/full"));
    assert(!matrix.save((work / "missing" / "matrix.bin").string()));

    StripePattern pattern;
    assert(pattern.parse("16"));
    std::filesystem::create_symlink("/dev/full", work / "failed_greedy.bin");
    assert(!generateSeedMatrix(pattern, PlacementStrategy::GREEDY, ErrorType::SUBSTITUTION,
                               (work / "failed").string(), 1, 61, 0, 0.999));
    assert(generateSeedMatrix(pattern, PlacementStrategy::GREEDY, ErrorType::SUBSTITUTION,
                              (work / "good").string(), 1, 61, 0, 0.999));
    assert(loaded.load((work / "good_greedy.bin").string()));
    assert(loaded.hasTruncMatrix());
    std::filesystem::remove_all(work);
    std::cout << "PASS: index boundaries and seed-matrix write failures\n";
}

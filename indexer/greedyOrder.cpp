#include "greedyOrder.hpp"

#include <cstdint>
#include <fstream>
#include <vector>

#include "SkippyConstants.hpp"
#include "stripe_pattern.h"

std::vector<int> buildGreedyOrder(int readLength, int span, const std::vector<int>& careOffsets,
                                  const std::vector<int>& protectOffsets);

std::vector<std::vector<uint32_t>> gGreedyOrders;

void initialiseGreedyOrders(const StripePattern& pattern)
{
    gGreedyOrders.resize(MAX_READ_LENGTH + 1);

    for (int L = 1; L <= MAX_READ_LENGTH; ++L) {
        std::vector<int> order =
            buildGreedyOrder(L, pattern.window, pattern.positions, pattern.protect_positions);

        gGreedyOrders[L] = std::vector<uint32_t>(order.begin(), order.end());
    }
}

bool loadGreedyOrders(const std::string& filename)
{
    std::ifstream in(filename, std::ios::binary);

    if (!in)
        return false;

    uint32_t version;
    uint32_t maxLength;

    in.read(reinterpret_cast<char*>(&version), sizeof(version));

    in.read(reinterpret_cast<char*>(&maxLength), sizeof(maxLength));

    if (version != 1)
        return false;

    gGreedyOrders.clear();
    gGreedyOrders.resize(maxLength);

    for (uint32_t i = 0; i < maxLength; ++i) {
        uint32_t count;

        in.read(reinterpret_cast<char*>(&count), sizeof(count));

        gGreedyOrders[i].resize(count);

        if (count) {
            in.read(reinterpret_cast<char*>(gGreedyOrders[i].data()), count * sizeof(uint32_t));
        }
    }

    return true;
}

bool saveGreedyOrders(const std::string& filename)
{
    std::ofstream out(filename, std::ios::binary);

    if (!out)
        return false;

    uint32_t version = 1;
    uint32_t maxLength = static_cast<uint32_t>(gGreedyOrders.size());

    out.write(reinterpret_cast<char*>(&version), sizeof(version));

    out.write(reinterpret_cast<char*>(&maxLength), sizeof(maxLength));

    for (const auto& vec : gGreedyOrders) {
        uint32_t count = static_cast<uint32_t>(vec.size());

        out.write(reinterpret_cast<char*>(&count), sizeof(count));

        if (count) {
            out.write(reinterpret_cast<const char*>(vec.data()), count * sizeof(uint32_t));
        }
    }

    return true;
}

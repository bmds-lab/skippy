#pragma once

#include <cstdint>
#include <string>
#include <vector>

struct StripePattern;

extern std::vector<std::vector<uint32_t>> gGreedyOrders;

void initialiseGreedyOrders(
    const StripePattern& pattern);

bool loadGreedyOrders(
    const std::string& filename);

bool saveGreedyOrders(
    const std::string& filename);

void initialiseGreedyOrders(
    const StripePattern& pattern);

bool loadGreedyOrders(
    const std::string& filename);

bool saveGreedyOrders(
    const std::string& filename);


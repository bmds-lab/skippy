// report_seed_matrix.hpp
#pragma once

#include "SeedMatrix.hpp"
#include "stripe_pattern.h"
#include "MonteCarlo.hpp"

#include <vector>
#include <string>

void reportSeedMatrix(
    SeedMatrix&                      sm,
    const std::vector<uint16_t>&     matrix,
    const std::vector<int>&          readLengths,
    const std::vector<double>&       errorRates,
    const StripePattern&             pattern,
    const std::string&               baseName,
    PlacementStrategy                strategy);

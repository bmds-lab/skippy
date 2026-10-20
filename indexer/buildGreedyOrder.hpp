#pragma once

#include <vector>

std::vector<int> buildGreedyOrder(int readLength, int span, const std::vector<int>& careOffsets,
                                  const std::vector<int>& protectOffsets);

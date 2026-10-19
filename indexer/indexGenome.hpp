/*
        INDEXGENOME.HPP
        ---------------
        indexKISS

        Created by Shlomo Geva on 13/7/2023.
        Updated: StripePattern replaces bool stripe.
*/
#pragma once

#include <math.h>
#include <string.h>
#include <sys/stat.h>

#include <map>
#include <vector>

#include "protected_vector.hpp"
#include "stripe_pattern.h"

char *read_entire_file(const char *filename, uint64_t &fileSize);
char *load_genome_file(const std::string &fastaFile,
                       std::map<uint32_t, std::string> &referenceIDMap, uint32_t &genomeSize);
void index_kmers(char *genome, uint64_t genomeSize,
                 std::vector<protected_vector<uint32_t>> &kmersMap,
                 uint32_t MASK,
                 const StripePattern &pat,
                 const StripePattern &pat2);

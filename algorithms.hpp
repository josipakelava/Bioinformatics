#pragma once

#include <cstdint>
#include <string>
#include <algorithm>
#include "kmer_util.hpp"
#include "lib/libbf/bf/all.hpp"


bool one_sided_kBF(kmer_t query, const bf::basic_bloom_filter& bf) {
    if(bf.lookup(query)) {
        vector<kmer_t> neighbours = {neighbor_left_set(query), neighbor_right_set(query)};
        return contains_set(neighbours, bf);
    }
    return false;
}


#pragma once

#include <cstdint>
#include <string>
#include <algorithm>
#include "kmer_util.hpp"
#include "lib/libbf/bf/all.hpp"


bool one_sided_kBF(kmer_t query, const bf::basic_bloom_filter& bf) {
    if(bf.lookup(query)) {
        vector<kmer_t> left(neighbor_left_set(query));
        vector<kmer_t> right(neighbor_right_set(query));
        return contains_set(left, bf) | contains_set(right, bf);
    }
    return false;
}


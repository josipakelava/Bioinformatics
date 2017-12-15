#pragma once

#include <cstdint>
#include <string>
#include <algorithm>
#include "kmer_util.hpp"
#include "lib/libbf/bf/all.hpp"
#include <unordered_set>


struct OneSided {
    bool one_sided_kBF(kmer_t query, const bf::basic_bloom_filter &bf) {
        if (bf.lookup(query)) {
            vector<kmer_t> left(neighbor_left_set(query));
            vector<kmer_t> right(neighbor_right_set(query));
            return contains_set(left, bf) | contains_set(right, bf);
        }
        return false;
    }

    bool operator()(kmer_t query, const bf::basic_bloom_filter& bf) {
        return one_sided_kBF(query, bf);
    }
};

struct TwoSided {

    unordered_set<kmer_t> edge_kmer;

    bool two_sided_kBF(kmer_t query, const bf::basic_bloom_filter& bf) {
        if (bf.lookup(query)) {
            vector<kmer_t> left(neighbor_left_set(query));
            vector<kmer_t> right(neighbor_right_set(query));

            bool contains_left = contains_set(left, bf);
            bool contains_right = contains_set(right, bf);

            if (contains_left && contains_right) return true;
            if (contains_left || contains_right)
                if(edge_kmer.find(query) != edge_kmer.end())
                    return true;
        }
        return false;
    }

    bool operator()(kmer_t query, const bf::basic_bloom_filter& bf) {
        return two_sided_kBF(query, bf);
    }
};


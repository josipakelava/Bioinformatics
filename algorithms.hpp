#pragma once

#include <cstdint>
#include <string>
#include <algorithm>
#include "kmer_util.hpp"
#include "lib/libbf/bf/all.hpp"
#include <unordered_set>

struct BasicBF {

    bf::basic_bloom_filter *bf;
    explicit BasicBF(const unordered_set<kmer_t> &set) {
        bf = new bf::basic_bloom_filter(bf::make_hasher(2), set.size()*10);
        for (kmer_t kmer : set) {
            bf->add(kmer);
        }
    }
    bool lookup(kmer_t query) {
        return bf->lookup(query) > 0;
    }

    ~BasicBF() {
        delete bf;
    }
};

struct OneSided {

    bf::basic_bloom_filter *bf;
    explicit OneSided(const unordered_set<kmer_t> &set) {

        bf = new bf::basic_bloom_filter(bf::make_hasher(2), set.size()*10);
        for (kmer_t kmer : set) {
            bf->add(kmer);
        }

    }
    bool lookup(kmer_t query) {
        if (bf->lookup(query)) {
            vector<kmer_t> left(neighbor_left_set(query));
            vector<kmer_t> right(neighbor_right_set(query));
            return contains_set(left, *bf) | contains_set(right, *bf);
        }
        return false;
    }

    ~OneSided() {
        delete bf;
    }
};

struct TwoSided {

    unordered_set<kmer_t> edge_kmer;

    bf::basic_bloom_filter *bf;
    TwoSided(const unordered_set<kmer_t> &set, const unordered_set<kmer_t> &edges) : edge_kmer(edges) {

        bf = new bf::basic_bloom_filter(bf::make_hasher(2), set.size()*10);
        for (kmer_t kmer : set) {
            bf->add(kmer);
        }
    }

    bool lookup(kmer_t query) {
        if (bf->lookup(query)) {
            vector<kmer_t> left(neighbor_left_set(query));
            vector<kmer_t> right(neighbor_right_set(query));

            bool contains_left = contains_set(left, *bf);
            bool contains_right = contains_set(right, *bf);

            if (contains_left && contains_right) return true;
            if (contains_left || contains_right)
                if(edge_kmer.find(query) != edge_kmer.end())
                    return true;
        }
        return false;
    }

    ~TwoSided() {
        delete bf;
    }
};

struct SparseKBF {

    bf::basic_bloom_filter *bf;
    int s;

    explicit SparseKBF(const unordered_set<kmer_t> &sparse_set, int s) : s(s) {

        bf = new bf::basic_bloom_filter(bf::make_hasher(2), sparse_set.size()*10);
        for (kmer_t kmer : sparse_set) {
            bf->add(kmer);
        }

    }

    bool lookup(kmer_t query) {
        if (bf->lookup(query)) {
            vector<kmer_t> left(neighbor_left_set(query));
            vector<kmer_t> right(neighbor_right_set(query));
            return contains_set(left, *bf) | contains_set(right, *bf);
        }
        return false;
    }

    bool lookupLeft(kmer_t query, int i) {

    }

    bool lookupRight(kmer_t query, int i) {

    }

    ~SparseKBF() {
        delete bf;
    }
};



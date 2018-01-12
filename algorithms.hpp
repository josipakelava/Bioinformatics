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
            return contains_set(neighbor_left_set(query), *bf) | contains_set(neighbor_right_set(query), *bf);
        }
        return false;
    }

    ~OneSided() {
        delete bf;
    }
};

bool decidePresent(kmer_t query, const unordered_set<kmer_t> &edges, bool containsLeft, bool containsRight) {
    if (containsRight && containsLeft) {
        return true;
    }

    if (containsRight || containsLeft) {
        if(edges.find(query) != edges.end()) {
            return true;
        }
    }

    return false;
}

struct TwoSided {

    unordered_set<kmer_t> edges;

    bf::basic_bloom_filter *bf;
    TwoSided(const unordered_set<kmer_t> &set, const unordered_set<kmer_t> &edges) : edges(edges) {

        bf = new bf::basic_bloom_filter(bf::make_hasher(2), set.size()*10);
        for (kmer_t kmer : set) {
            bf->add(kmer);
        }
    }

    bool lookup(kmer_t query) {
        if (bf->lookup(query)) {
            return decidePresent(query, edges, contains_set(neighbor_left_set(query), *bf), contains_set(neighbor_right_set(query), *bf));
        }

        return false;
    }

    ~TwoSided() {
        delete bf;
    }
};

struct SparseKBF {

    bf::basic_bloom_filter *bf;
    unordered_set<kmer_t> edges;
    int s;

    explicit SparseKBF(const unordered_set<kmer_t> &sparse_set, const unordered_set<kmer_t> &edges, int s) : edges(edges), s(s) {

        bf = new bf::basic_bloom_filter(bf::make_hasher(2), sparse_set.size()*10);
        for (kmer_t kmer : sparse_set) {
            bf->add(kmer);
        }

    }

    bool lookup(kmer_t query) {
        if (bf->lookup(query)) {
            if (strictContainsNeighbours(query, s, s)) {
                return true;
            }
        }

        for(int i = 0; i < s - 1; i++) {
            if (strictContainsNeighbours(query, i, s - (i + 1))) {
                return true;
            }
        }

        return false;
    }

    bool strictContainsNeighbours(kmer_t query, int left, int right) {
        return decidePresent(query, edges, contains_set(strict_neighbor_set(query, left, 0), *bf), contains_set(
                strict_neighbor_set(query, 0, right), *bf));
    }

    ~SparseKBF() {
        delete bf;
    }
};

struct SparseRelaxedKBF {

    bf::basic_bloom_filter *bf;
    unordered_set<kmer_t> edges;
    int s;

    explicit SparseRelaxedKBF(const unordered_set<kmer_t> &sparse_set, const unordered_set<kmer_t> &edges, int s) : edges(edges), s(s) {

        bf = new bf::basic_bloom_filter(bf::make_hasher(2), sparse_set.size()*10);
        for (kmer_t kmer : sparse_set) {
            bf->add(kmer);
        }

    }

    bool lookup(kmer_t query) {
        if (bf->lookup(query)) {
            if (relaxedContainsNeighbours(query, s+1, s+1)) {
                return true;
            }
        } else {
            for(int i = 0; i < s; i++) {
                if (relaxedContainsNeighbours(query, i+1, s - i)) {
                    return true;
                }
            }
        }
        return false;
    }

    bool relaxedContainsNeighbours(kmer_t query, int left, int right) {
        return decidePresent(query, edges, contains_set(relaxed_neighbor_set(query, left, 0), *bf), contains_set(relaxed_neighbor_set(query, 0, right), *bf));
    }


    ~SparseRelaxedKBF() {
        delete bf;
    }
};




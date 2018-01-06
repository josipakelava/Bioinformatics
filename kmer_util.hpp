#ifndef BIOINFORMATIKA_PROJEKT_KMER_UTIL_HPP
#define BIOINFORMATIKA_PROJEKT_KMER_UTIL_HPP

#include <cstdint>
#include <string>
#include <algorithm>
#include <unordered_set>
#include <queue>
#include "DEBUG_UTIL.hpp"
#include "lib/libbf/bf/all.hpp"

using namespace std;

const int KMER_LENGTH = 20; //isto kao i u radu
const int KMER_SHIFT_LEFT = (2 * (KMER_LENGTH - 1));
using kmer_t = uint64_t;

const kmer_t INVALID_KMER = UINT64_MAX;
const kmer_t KMER_MASK = ((uint64_t) 1 << (2 * KMER_LENGTH)) - 1;

const int QUERY_SET_SIZE = 1e4;

inline kmer_t base_to_bits(char c) {
    if (c == 'A' || c == 'a') return 0;
    if (c == 'T' || c == 't') return 1;
    if (c == 'G' || c == 'g') return 2;
    if (c == 'C' || c == 'c') return 3;

    return INVALID_KMER;
}

inline char bits_to_base(kmer_t kmer) {
    const char m[4] = {'A', 'T', 'G', 'C'};
    return m[kmer & 3];
}

kmer_t substring_to_kmer(const string &s, int pos) {
    if (pos + KMER_LENGTH > s.length()) {
        return INVALID_KMER;
    }

    kmer_t kmer = 0;
    for (int i = 0; i < KMER_LENGTH; i++) {
        kmer = kmer << 2 | base_to_bits(s[i + pos]);
    }

    return kmer;
}

kmer_t string_to_kmer(const string &s) {
    return substring_to_kmer(s, 0);
}

string kmer_to_string(kmer_t kmer) {
    string s = "";
    for (int i = 0; i < KMER_LENGTH; i++) {
        s += bits_to_base(kmer);
        kmer >>= 2;
    }
    reverse(s.begin(), s.end());
    return s;
}

inline kmer_t add_base_right(kmer_t kmer, char base) {
    return ((kmer << 2) & KMER_MASK) | base_to_bits(base);
}

inline kmer_t add_base_left(kmer_t kmer, char base) {
    return (kmer >> 2) | base_to_bits(base) << KMER_SHIFT_LEFT;
}

bool contains_set(const vector<kmer_t> &query, const bf::basic_bloom_filter &bf) {
    for (kmer_t kmer : query) {
        if (bf.lookup(kmer)) {
            return true;
        }
    }
    return false;
}

unordered_set<kmer_t> query_set(const vector<kmer_t> &kmers) {
    unordered_set<kmer_t> query_set;

    mt19937 random;
    random.seed(std::random_device()());

    uniform_int_distribution<> kmers_dist(0, kmers.size() - 1);
    uniform_int_distribution<> base_dist(0, 2*KMER_LENGTH - 1);

    auto it = kmers.begin();
    for(int i = 0; i < QUERY_SET_SIZE; i++) {
        kmer_t original = *(it + kmers_dist(random));
        kmer_t mutated = original ^ 1 << base_dist(random);
        query_set.insert(mutated);
    }

    return query_set;
}

unordered_set<kmer_t> generate_set(const vector<string> &sequnces){

    unordered_set<kmer_t> set;
    for (string seq : sequnces) {
        kmer_t kmer = string_to_kmer(seq);
        set.insert(kmer);
        for (int i = KMER_LENGTH; i < seq.length(); i++) {
            kmer = add_base_right(kmer, seq[i]);
            set.insert(kmer);
        }
    }
    return set;
}

unordered_set<kmer_t> generate_set_edge(const vector<string> &sequnces, unordered_set<kmer_t> &edge_kmer){

    unordered_set<kmer_t> set;
    for (string seq : sequnces) {
        kmer_t kmer = string_to_kmer(seq);
        set.insert(kmer);
        edge_kmer.insert(kmer);
        for (int i = KMER_LENGTH; i < seq.length(); i++) {
            kmer = add_base_right(kmer, seq[i]);
            set.insert(kmer);
        }
        edge_kmer.insert(kmer);
    }
    return set;
}

unordered_set<kmer_t> generate_sparse_set(const string &seq, const int s){

    unordered_set<kmer_t> set;
    kmer_t kmer;

    for (int i = 0; i < seq.length(); i+= s + 1) {
        kmer = substring_to_kmer(seq, i);
        set.insert(kmer);
    }

    return set;
}

vector<kmer_t> neighbor_left_set(kmer_t query) {
    vector<kmer_t> neighbors = {add_base_left(query, 'A'), add_base_left(query, 'G'), add_base_left(query, 'T'),
                                add_base_left(query, 'C')};
    return neighbors;
}

vector<kmer_t> neighbor_right_set(kmer_t query) {
    vector<kmer_t> neighbors = {add_base_right(query, 'A'), add_base_right(query, 'G'), add_base_right(query, 'T'),
                                add_base_right(query, 'C')};
    return neighbors;
}

vector<kmer_t> neighbor_set(kmer_t query, int left, int right) {
    vector<kmer_t> neighbors;

    queue<kmer_t> result;
    result.push(query);

    queue<kmer_t> current;
    for (int i = 0; i < left; ++i) {

        while(result.size() > 0) {
            kmer_t kmer = result.front(); result.pop();
            current.push(kmer);
        }

        while(current.size() > 0) {
            kmer_t kmer = current.front(); current.pop();
            vector<kmer_t> neighbors = neighbor_left_set(kmer);
            for(kmer_t neighbor : neighbors) {
                result.push(neighbor);
            }
        }
    }

    while(result.size() > 0) {
        kmer_t kmer = result.front(); result.pop();
        neighbors.push_back(kmer);
    }

    result.push(query);

    for (int i = 0; i < right; ++i) {

        while(result.size() > 0) {
            kmer_t kmer = result.front(); result.pop();
            current.push(kmer);
        }

        while(current.size() > 0) {
            kmer_t kmer = current.front(); current.pop();
            vector<kmer_t> neighbors = neighbor_right_set(kmer);
            for(kmer_t neighbor : neighbors) {
                result.push(neighbor);
            }
        }
    }

    while(result.size() > 0) {
        kmer_t kmer = result.front(); result.pop();
        neighbors.push_back(kmer);
    }

    return neighbors;
}

vector<kmer_t> neighbor_right_set(kmer_t query, const int s) {
    vector<kmer_t> neighbors = {add_base_right(query, 'A'), add_base_right(query, 'G'), add_base_right(query, 'T'),
                                add_base_right(query, 'C')};
    return neighbors;
}

#endif //BIOINFORMATIKA_PROJEKT_KMER_UTIL_HPP

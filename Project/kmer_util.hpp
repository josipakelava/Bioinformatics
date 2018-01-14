    #ifndef BIOINFORMATIKA_PROJEKT_KMER_UTIL_HPP
#define BIOINFORMATIKA_PROJEKT_KMER_UTIL_HPP

#include <cstdint>
#include <string>
#include <algorithm>
#include <unordered_set>
#include <queue>
#include <vector>
#include "DEBUG_UTIL.hpp"
#include "../lib/libbf/bf/all.hpp"
#include <unordered_map>
#include <set>
#include <iostream>
#include <assert.h>

using namespace std;

using kmer_t = uint64_t;
const int KMER_LENGTH = 20; //isto kao i u radu
const kmer_t KMER_SHIFT_LEFT = (uint64_t)2 * KMER_LENGTH - 2;

const kmer_t INVALID_KMER = UINT64_MAX;
const kmer_t KMER_MASK = ((uint64_t) 1 << (2 * KMER_LENGTH)) - 1;

const int QUERY_SET_SIZE = 1e4;

/**
 * Converts base from char to kmer
 * @param c base
 * @return kmer for base
 */
inline kmer_t base_to_bits(char c) {
    if (c == 'A' || c == 'a') return 0;
    if (c == 'T' || c == 't') return 1;
    if (c == 'G' || c == 'g') return 2;
    if (c == 'C' || c == 'c') return 3;

    assert(true);
    return INVALID_KMER;
}

/**
* Converts base from kmer to char
* @param kmer base
* @return bits for base
*/
inline char bits_to_base(kmer_t kmer) {
    const char m[4] = {'A', 'T', 'G', 'C'};
    return m[kmer & 3];
}

/**
 * Converts part of string to  kmer
 * @param s string to convert
 * @param pos start of convertions
 * @return converted string
 */
kmer_t substring_to_kmer(const string &s, int pos) {
    if (pos + KMER_LENGTH > s.length()) {
        cout << "INVALID" << endl;
        return INVALID_KMER;
    }

    kmer_t kmer = 0;
    for (int i = 0; i < KMER_LENGTH; i++) {
        kmer = kmer << 2 | base_to_bits(s[i + pos]);
    }

    if ((kmer&KMER_MASK) != kmer) cout << "ERROR!" << endl;

    return kmer&KMER_MASK;
}

/**
 * Converts string to  kmer
 * @param s string to convert
 * @return converted string
 */
kmer_t string_to_kmer(const string &s) {
    return substring_to_kmer(s, 0);
}

/**
 * Converts kmer to string
 * @param kmer kmer to convert
 * @return converted kmer
 */
string kmer_to_string(kmer_t kmer) {
    string s = "";
    for (int i = 0; i < KMER_LENGTH; i++) {
        s += bits_to_base(kmer);
        kmer >>= 2;
    }
    reverse(s.begin(), s.end());
    return s;
}

/**
 * Adds base to the right side of kmer
 * @param kmer
 * @param base base to add to the kmer
 * @return new kmer
 */
inline kmer_t add_base_right(kmer_t kmer, char base) {
    return ((kmer << 2) & KMER_MASK) | base_to_bits(base);
}

/**
* Adds base to the left side of kmer
* @param kmer
* @param base base to add to the kmer
* @return new kmer
*/
inline kmer_t add_base_left(kmer_t kmer, char base) {
    return (kmer >> 2) | (base_to_bits(base) << KMER_SHIFT_LEFT);
}

/**
 * Checks if bloom filter contains at least one kmer from query set
 * @param query set of kmers
 * @param bf bloom filter
 * @return true if at least one kmer exist
 */
bool contains_set(const unordered_set<kmer_t> &query, const bf::basic_bloom_filter &bf) {
    for (kmer_t kmer : query) {
        if (bf.lookup(kmer)) {
            return true;
        }
    }
    return false;
}

/**
 * Generates set of queries. Number of queries is defined by parameter QUERY_SET_SIZE.
 * @param kmers set of possible kmers for query
 * @return set of queries
 */
unordered_set<kmer_t> query_set(const unordered_set<kmer_t> &kmers) {
    unordered_set<kmer_t> query_set;

    mt19937 random;
    random.seed(std::random_device()());

    uniform_int_distribution<> base_dist(0, 2 * KMER_LENGTH - 1);

    auto it = kmers.begin();
    for(int i = 0; i < QUERY_SET_SIZE; i++) {
        kmer_t original = *it;
        kmer_t mutated = original ^ 1 << base_dist(random);
        query_set.insert(mutated);
        it++;
    }

    return query_set;
}

/**
 * Generates kmers from sequnecs
 * @param sequnces sequnces to convert
 * @return set of kmers
 */
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

/**
 * Generates kmers from sequnces while saving edge kmers
 * @param sequnces sequnces to convert
 * @param edge_kmer refernce to set of edge kmers
 * @param numberOfPotentialEdges 
 * @return set of generated kmers
 */
unordered_set<kmer_t> generate_set_with_edges(const vector<string> &sequnces, unordered_set<kmer_t> &edge_kmer, int &numberOfPotentialEdges){

    unordered_set<kmer_t> kmers;
    unordered_set<kmer_t> potential_edge_kmers;

    for (string seq : sequnces) {
        kmer_t kmer = string_to_kmer(seq);
        potential_edge_kmers.insert(kmer);
        int i = KMER_LENGTH;
        for (; i < seq.length() - 1; i++) {
            kmer = add_base_right(kmer, seq[i]);
            kmers.insert(kmer);
        }
        kmer = add_base_right(kmer, seq[i]);
        potential_edge_kmers.insert(kmer);
    }

    numberOfPotentialEdges = potential_edge_kmers.size();

    for(kmer_t kmer : potential_edge_kmers) {
        if (kmers.find(kmer) == kmers.end()) {
            edge_kmer.insert(kmer);
        }
    }

    kmers.insert(edge_kmer.begin(), edge_kmer.end());

    return kmers;
}

unordered_set<kmer_t> generate_sparse_set(const string &seq, unordered_set<kmer_t> &edge_kmer, const int s){

    unordered_set<kmer_t> set;
    kmer_t kmer;

    edge_kmer.insert(string_to_kmer(seq));

    int i = 0;
    for (; i < seq.length() - KMER_LENGTH + 1; i += s + 1) {
        kmer = substring_to_kmer(seq, i);
        set.insert(kmer);
    }
    for(i -= s + 1; i < seq.length() - KMER_LENGTH+1; i++)
        edge_kmer.insert(substring_to_kmer(seq, i));
    return set;
}

unordered_set<kmer_t> generate_best_fit_set(const vector<string> &sequences, unordered_set<kmer_t> &edge_kmer, int s) {
    unordered_set<kmer_t> kmers;
    const int MOD = s + 1;
    for(auto& seq : sequences) {
        int count[MOD];
        int mx_ind = 0;
        for(int i = 0; i < MOD; i++) count[i] = 0;
        for(int i = 0; i < seq.length() - KMER_LENGTH + 1; i++) {
            if(kmers.find(substring_to_kmer(seq, i)) != kmers.end()) {
                const int ind = i % MOD;
                if(++count[ind] > count[mx_ind])
                    mx_ind = ind;
            }
        }

        int i = 0;
        for(; i <= mx_ind; i++) edge_kmer.insert(substring_to_kmer(seq, i));
        for (i = mx_ind; i < seq.length() - KMER_LENGTH + 1; i += MOD) {
            kmers.insert(substring_to_kmer(seq, i));
        }
        for(i -= MOD; i < seq.length() - KMER_LENGTH + 1; i++)
            edge_kmer.insert(substring_to_kmer(seq, i));
    }
    return kmers;
}

struct KmerNode {
    kmer_t kmer;
    int deg;
    unordered_set<unordered_set<kmer_t>*> neighbours;

    KmerNode(kmer_t kmer) : kmer(kmer), deg(0) {}

    void add(unordered_set<kmer_t>* n) {
        neighbours.insert(n);
        deg++;
    }

    void remove(unordered_set<kmer_t>* n) {
        if(deg == 0) {
            cout << "deg < 0" << endl;
            return;
        }
        neighbours.erase(n);
        deg--;
    }

};

struct NodeCMP {
    bool operator()(const KmerNode* first, const KmerNode* second) {
        if (first->deg == second->deg) {
            return first->kmer > second->kmer;
        }

        return first->deg > second->deg;
    }

};

void compute_graph(set<KmerNode*, NodeCMP>& s, unordered_map<kmer_t, KmerNode*>& kmerToPtr,
    vector<KmerNode*>& allPtrs, unordered_map<kmer_t, unordered_set<kmer_t>>& m) {
    for(auto& kv : m) {
        for(auto kmer : kv.second) {
            if(kmerToPtr.find(kmer) == kmerToPtr.end()) {
                KmerNode *node = new KmerNode(kmer);
                kmerToPtr[kmer] = node;
                node->add(&kv.second);
                s.insert(node);
                allPtrs.push_back(node);
            } else {
                s.erase(kmerToPtr[kmer]);
                kmerToPtr[kmer]->add(&kv.second);
                s.insert(kmerToPtr[kmer]);
            }
        }
    }
}

unordered_set<kmer_t> generate_hitting_set_kmers(const vector<string>& sequences, unordered_set<kmer_t> &edge_kmer) {
    unordered_map<kmer_t, unordered_set<kmer_t>> leftN, rightN;
    for(const auto &seq : sequences) {
        kmer_t left = string_to_kmer(seq);
        kmer_t mid = add_base_right(left, seq[KMER_LENGTH]);
        kmer_t right = add_base_right(mid, seq[KMER_LENGTH+1]);
        edge_kmer.insert(left); //first
        rightN[left].insert(mid);
        leftN[left].insert(left);
        rightN[left].insert(left);
        for(int i = KMER_LENGTH+2; i < seq.size(); i++) {
            leftN[mid].insert(left);
            rightN[mid].insert(right);
            leftN[mid].insert(mid);
            rightN[mid].insert(mid);
            left = mid;
            mid = right;
            right = add_base_right(right, seq[i]);
        }
        leftN[right].insert(right);
        rightN[right].insert(right);
        leftN[right].insert(mid);
        edge_kmer.insert(right);
    }

    set<KmerNode*, NodeCMP> s;
    unordered_map<kmer_t, KmerNode*> kmerToPtr;
    vector<KmerNode*> allPtrs;

    compute_graph(s, kmerToPtr, allPtrs, leftN);
    compute_graph(s, kmerToPtr, allPtrs, rightN);

    unordered_set<kmer_t> kmers;
    while(!s.empty()) {
        KmerNode* mx = *s.begin();
        s.erase(mx);
        kmers.insert(mx->kmer);
        for(auto ns : mx->neighbours) {
            for(auto kmer : *ns) {
                if(kmer == mx->kmer || kmerToPtr.find(kmer) == kmerToPtr.end()) continue;
                s.erase(kmerToPtr[kmer]);
                kmerToPtr[kmer]->remove(ns);
                if (kmerToPtr[kmer]->deg > 0) {
                    s.insert(kmerToPtr[kmer]);
                } else {
                    kmerToPtr.erase(kmer);
                }
            }
        }
    }

    while(allPtrs.size()) {
        auto p = allPtrs.back();
        allPtrs.pop_back();
        delete p;
    }

    return kmers;
}

/**
 * Generates left neighbours for query
 * @param query
 * @return set of left neighbours
 */
unordered_set<kmer_t> neighbor_left_set(kmer_t query) {
    unordered_set<kmer_t> neighbors = {add_base_left(query, 'A'), add_base_left(query, 'G'), add_base_left(query, 'T'),
                                add_base_left(query, 'C')};
    return neighbors;
}

/**
* Generates right neighbours for query
* @param query
* @return set of right neighbours
*/
unordered_set<kmer_t> neighbor_right_set(kmer_t query) {
    unordered_set<kmer_t> neighbors = {add_base_right(query, 'A'), add_base_right(query, 'G'), add_base_right(query, 'T'),
                                add_base_right(query, 'C')};
    return neighbors;
}

/**
* Generates neighbour set for strict k-mer sparsification
* @param query
* @param left number of left k-mers
* @param right number of right k-mers
* @return set of strict neighbours
 */
unordered_set<kmer_t> strict_neighbor_set(kmer_t query, int left, int right) {
    unordered_set<kmer_t> neighbors;

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
            unordered_set<kmer_t> neighbors = neighbor_left_set(kmer);
            for(kmer_t neighbor : neighbors) {
                result.push(neighbor);
            }
        }
    }

    while(result.size() > 0 && left > 0) {
        kmer_t kmer = result.front(); result.pop();
        neighbors.insert(kmer);
    }

    result.push(query);

    for (int i = 0; i < right; ++i) {

        while(result.size() > 0) {
            kmer_t kmer = result.front(); result.pop();
            current.push(kmer);
        }

        while(current.size() > 0) {
            kmer_t kmer = current.front(); current.pop();
            unordered_set<kmer_t> neighbors = neighbor_right_set(kmer);
            for(kmer_t neighbor : neighbors) {
                result.push(neighbor);
            }
        }
    }

    while(result.size() > 0 && right > 0) {
        kmer_t kmer = result.front(); result.pop();
        neighbors.insert(kmer);
    }

    return neighbors;
}

/**
 * Generates neighbour set for relaxed k-mer sparsification
 * @param query
 * @param left number of left k-mers
 * @param right number of right k-mers
 * @return set of relaxed neighbours
 */
unordered_set<kmer_t> relaxed_neighbor_set(kmer_t query, int left, int right) {
    unordered_set<kmer_t> neighborsFinal;

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
            unordered_set<kmer_t> neighbors = neighbor_left_set(kmer);
            for(kmer_t neighbor : neighbors) {
                result.push(neighbor);
                neighborsFinal.insert(neighbor);
            }
        }
    }

    queue<kmer_t >().swap(result);
    result.push(query);

    for (int i = 0; i < right; ++i) {

        while(result.size() > 0) {
            kmer_t kmer = result.front(); result.pop();
            current.push(kmer);
        }

        while(current.size() > 0) {
            kmer_t kmer = current.front(); current.pop();
            unordered_set<kmer_t> neighbors = neighbor_right_set(kmer);
            for(kmer_t neighbor : neighbors) {
                result.push(neighbor);
                neighborsFinal.insert(neighbor);
            }
        }
    }

    return neighborsFinal;
}

#endif //BIOINFORMATIKA_PROJEKT_KMER_UTIL_HPP

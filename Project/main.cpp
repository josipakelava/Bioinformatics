#include <iostream>
#include "../lib/libbf/bf/all.hpp"
#include "fasta_io.hpp"
#include "kmer_util.hpp"
#include "algorithms.hpp"

using namespace std;

void kmer_test(int argc, char **argv) {
    if (argc < 2) return;

    vector<string> reads;
    read_fasta(argv[1], reads);
    int i = 0;
    for(const auto& s : reads) {
        cout << "Read " << i++ << endl;
        cout << s << endl;
    }

    kmer_t kmer = string_to_kmer(reads[0]);
    cout << "kmer1: " << kmer_to_string(kmer) << endl;
    cout << "kmer2: " << kmer_to_string(add_base_right(kmer, 'A')) << endl;
    cout << "kmer2: " << kmer_to_string(add_base_left(kmer, 'A')) << endl;
}

void bloom_test() {
    bf::basic_bloom_filter b(0.8, 100);
    // Add two elements.
    b.add("foo");
    b.add(42);

    // Test set membership
    std::cout << b.lookup("foo") << std::endl;  // 1
    std::cout << b.lookup("bar") << std::endl;  // 0
    std::cout << b.lookup(42) << std::endl;     // 1

    // Remove all elements.
    b.clear();
    std::cout << b.lookup("foo") << std::endl;  // 0
    std::cout << b.lookup(42) << std::endl;     // 0
}

void contains_set_test() {
    cout << "Contains set" << std::endl;
    unordered_set<kmer_t> query = {1, 2, 3, 4, 5};
    bf::basic_bloom_filter bf(0.8, 100);
    cout << contains_set(query, bf) << std::endl; // 0
    bf.add(3);
    bf.add(5);
    cout << contains_set(query, bf) << std::endl; // 1
}

void test1(const unordered_set<kmer_t> &set) {
    BasicBF basicBF(set);
    OneSided kbf1(set);

    auto it = set.begin();
    it++;
    kmer_t random2 = *(it);
    kmer_t random1 = random2 ^(15 << 5);
    cout << "Test kbf1" << endl;
    cout << kbf1.lookup(random2) << endl;
    cout << kbf1.lookup(random1) << endl;

    cout << basicBF.lookup(random2) << endl;
    cout << basicBF.lookup(random1) << endl;

    cout << (set.find(random2) != set.end()) << endl;
    cout << (set.find(random1) != set.end()) << endl;
}

void test2(const unordered_set<kmer_t> &set, const unordered_set<kmer_t> &edges) {

    BasicBF basicBF(set);
    TwoSided kbf2(set, edges);

    auto it = set.begin();
    it++;
    kmer_t random2 = *(it);
    kmer_t random1 = random2 ^(15 << 5);
    cout << "Test kbf2" << endl;
    cout << kbf2.lookup(random2) << endl;
    cout << kbf2.lookup(random1) << endl;

    cout << basicBF.lookup(random2) << endl;
    cout << basicBF.lookup(random1) << endl;

    cout << (set.find(random2) != set.end()) << endl;
    cout << (set.find(random1) != set.end()) << endl;
}

void test3(const unordered_set<kmer_t> &set, const unordered_set<kmer_t> &edges, int s) {

    BasicBF basicBF(set);
    SparseRelaxedKBF kbf(set, edges, s);

    auto it = set.begin();
    it++;
    kmer_t random2 = *(it);
    kmer_t random1 = random2 ^(4 << 5);
    cout << "Test sparse" << endl;
    cout << kmer_to_string(random2) << endl;
    cout << kbf.lookup(random2) << endl;
    cout << kbf.lookup(random1) << endl;

    cout << basicBF.lookup(random2) << endl;
    cout << basicBF.lookup(random1) << endl;

    cout << (set.find(random2) != set.end()) << endl;
    cout << (set.find(random1) != set.end()) << endl;
}

int main(int argc, char **argv) {

//    vector<string> sequences;
//    read_fasta(argv[1], sequences);
//    unordered_set<kmer_t> kmers(generate_set(sequences));
//
//    unordered_set<kmer_t> edge_kmers;
//    unordered_set<kmer_t> kmers2(generate_set_edge(sequences, edge_kmers));
//
////    kmer_test(argc, argv);
//    test1(kmers);
//    test2(kmers, edge_kmers);
//
//    contains_set_test();

//    unordered_set<kmer_t> set = strict_neighbor_set(string_to_kmer("AAAAAAAAAAAAAAAAAAAA"), 0, 3);
//    for(kmer_t neighbour : set) {
//        cout << kmer_to_string(neighbour) << endl;
//    }
//    cout << set.size()<< endl;

//    unordered_set<kmer_t> set = relaxed_neighbor_set(string_to_kmer("AAAAAAAAAAAAAAAAAATA"), 0, 3);
//    for(kmer_t neighbour : set) {
//        cout << kmer_to_string(neighbour) << endl;
//    }
//    cout << set.size()<< endl;

    vector<string> sequences;
    read_fasta(argv[1], sequences);
    cout << "seq: " << sequences.size() << endl;
    unordered_set<kmer_t> edge_kmers;
    unordered_set<kmer_t> kmers(generate_hitting_set_kmers(sequences, edge_kmers));
    test3(kmers, edge_kmers, 1);

    return 0;
}
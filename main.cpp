#include <iostream>
#include "lib/libbf/bf/all.hpp"
#include "fasta_io.hpp"
#include "kmer_util.hpp"

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
    vector<kmer_t> query = {1, 2, 3, 4, 5};
    bf::basic_bloom_filter bf(0.8, 100);
    cout << contains_set(query, bf) << std::endl; // 0
    bf.add(3);
    bf.add(5);
    cout << contains_set(query, bf) << std::endl; // 1
}

int main(int argc, char **argv) {
//    kmer_test(argc, argv);
//    bloom_test();
    contains_set_test();
    return 0;
}
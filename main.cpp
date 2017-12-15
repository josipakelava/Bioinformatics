#include <iostream>
#include "lib/libbf/bf/all.hpp"
#include "fasta_io.hpp"
#include "kmer_util.hpp"
#include "algorithms.hpp"
#include <unordered_set>

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
void test1(const char* filePath) {
    vector<string> sequnces;
    read_fasta(filePath, sequnces);
    unordered_set<kmer_t> set;
    for (string seq : sequnces) {
        kmer_t kmer = string_to_kmer(seq);
        set.insert(kmer);
        for (int i = KMER_LENGTH; i < seq.length(); i++) {
            kmer = add_base_right(kmer, seq[i]);
            set.insert(kmer);
        }
    }
    bf::basic_bloom_filter bf(bf::make_hasher(2), set.size()*10);
    for (kmer_t kmer : set) {
        bf.add(kmer);
    }

    auto it = set.begin();
    it++;
    kmer_t random2 = *(it);
    kmer_t random1 = random2 ^(15 << 5);
    cout << "Test kbf1" << endl;
    cout << one_sided_kBF(random2, bf) << endl;
    cout << one_sided_kBF(random1, bf) << endl;

    cout << bf.lookup(random2) << endl;
    cout << bf.lookup(random1) << endl;

    cout << (set.find(random2) != set.end()) << endl;
    cout << (set.find(random1) != set.end()) << endl;
}

int main(int argc, char **argv) {
//    kmer_test(argc, argv);
    test1(argv[1]);
    contains_set_test();
    return 0;
}
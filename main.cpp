#include <iostream>

#include "fasta_io.h"
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
    cout << "kmer: " << kmer_to_string(kmer) << endl;
}

int main(int argc, char **argv) {
    kmer_test(argc, argv);
    return 0;
}
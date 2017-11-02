#include <iostream>

#include "fasta_io.h"

using namespace std;

void fasta_test(int argc, char **argv) {
    if (argc < 2) return;

    vector<string> reads;
    read_fasta(argv[1], reads);
    int i = 0;
    for(const auto& s : reads) {
        cout << "Read " << i++ << endl;
        cout << s << endl;
    }
}

int main(int argc, char **argv) {
    fasta_test(argc, argv);
    return 0;
}
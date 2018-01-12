#include <iostream>
#include "lib/libbf/bf/all.hpp"
#include "fasta_io.hpp"
#include "kmer_util.hpp"
#include "algorithms.hpp"

using namespace std;


void testFPR(const vector<string> sequences) {

    unordered_set<kmer_t> allKmersEdges;
    unordered_set<kmer_t> allKmers(generate_set_with_edges(sequences, allKmersEdges));

    BasicBF basicBF(allKmers);
    OneSided oneSidedBF(allKmers);
    TwoSided twoSidedBF(allKmers, allKmersEdges);

    unordered_set<kmer_t> bestFitKmersEdges;
    unordered_set<kmer_t> bestFitKmers(generate_best_fit_set(sequences, bestFitKmersEdges, 1));
    SparseKBF bestFitBF(bestFitKmers, bestFitKmersEdges, 1);


    unordered_set<kmer_t> hittingSetEdges;
    unordered_set<kmer_t> hittingSetKmers(generate_hitting_set_kmers(sequences, hittingSetEdges));
    SparseRelaxedKBF hittingSetBF(hittingSetKmers, hittingSetEdges, 1);

    unordered_set<kmer_t> querySet = query_set(allKmers);

    bool real = false;

    double falsePositiveNumBasic = 0;
    double falsePositiveNumOneSided = 0;
    double falsePositiveNumTwoSided = 0;
    double falsePositiveNumBestFit = 0;
    double falsePositiveNumHittingSet = 0;

    for (kmer_t kmer : querySet) {
        real = allKmers.find(kmer) != allKmers.end();

        if (!real && basicBF.lookup(kmer)) {
            falsePositiveNumBasic++;
        }

        if (!real && oneSidedBF.lookup(kmer)) {
            falsePositiveNumOneSided++;
        }

        if (!real && twoSidedBF.lookup(kmer)) {
            falsePositiveNumTwoSided++;
        }

        if (!real && bestFitBF.lookup(kmer)) {
            falsePositiveNumBestFit++;
        }

        if (!real && hittingSetBF.lookup(kmer)) {
            falsePositiveNumHittingSet++;
        }

//        if (real && !bestFitBF.lookup(kmer)) {
//            cout << "BF BF " << endl;
//        }
//
//        if (real && !hittingSetBF.lookup(kmer)) {
//            cout << "HS FN " << endl;
//        }
    }

    cout << "BASIC FPR " << falsePositiveNumBasic / querySet.size() << endl;
    cout << "ONE SIDED FPR " << falsePositiveNumOneSided / querySet.size() << endl;
    cout << "TWO SIDED FPR " << falsePositiveNumTwoSided / querySet.size() << endl;
    cout << "BEST FIT FPR " << falsePositiveNumBestFit/ querySet.size() << endl;
    cout << "HITTING SET FPR " << falsePositiveNumHittingSet / querySet.size() << endl;
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

    vector<string> sequences;
    read_fasta(argv[1], sequences);

    testFPR(sequences);

    return 0;
}
#include <iostream>
#include <ctime>
#include <chrono>
#include "../lib/libbf/bf/all.hpp"
#include "fasta_io.hpp"
#include "kmer_util.hpp"
#include "algorithms.hpp"

using namespace std;


void testFPR(const vector<string> sequences) {

    auto basicStart = chrono::system_clock::now();
    unordered_set<kmer_t> allKmers(generate_set(sequences));
    BasicBF basicBF(allKmers);
    auto basicEnd = chrono::system_clock::now();

    OneSided oneSidedBF(allKmers);

    auto twoSidedStart = chrono::system_clock::now();;
    unordered_set<kmer_t> allKmersEdges;
    int numberOfPotentialEdges = 0;
    allKmers = generate_set_with_edges(sequences, allKmersEdges, numberOfPotentialEdges);
    TwoSided twoSidedBF(allKmers, allKmersEdges);
    auto twoSidedEnd = chrono::system_clock::now();;

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

    cout << "Table 2. False Positive Rates" << endl;
    cout << "BASIC FPR " << falsePositiveNumBasic / querySet.size() << endl;
    cout << "ONE SIDED FPR " << falsePositiveNumOneSided / querySet.size() << endl;
    cout << "TWO SIDED FPR " << falsePositiveNumTwoSided / querySet.size() << endl;
    cout << "BEST FIT FPR " << falsePositiveNumBestFit/ querySet.size() << endl;
    cout << "HITTING SET FPR " << falsePositiveNumHittingSet / querySet.size() << endl;

    cout << endl;

    cout << "Table 3. Two-Sided k-mer Bloom Filter Overhead" << endl;
    cout << "No. of k-mers " << allKmers.size() << endl;
    cout << "No. of edge k-mers " << allKmersEdges.size() << endl;
    cout << "No. of potential edge k-mers " <<  numberOfPotentialEdges << endl;
//    nije dobro
    cout << "Initialization time (fold change) " << log2((twoSidedEnd - twoSidedStart).count()) - log2((basicEnd - basicStart).count())  << endl;

    cout << endl;

    cout << "Table 4. Number of k-mers Selected by Sparse k-mer Bloom Filter" << endl;
    cout << "No. of k-mers classic " << allKmers.size() << endl;
    cout << "No. of k-mers best match " << bestFitKmers.size() << endl;
//    nije dobro
    cout << "No. of k-mers hitting set " <<  hittingSetKmers.size() << endl;
}


int main(int argc, char **argv) {

    vector<string> sequences;
    read_fasta(argv[1], sequences);

    testFPR(sequences);

    return 0;
}
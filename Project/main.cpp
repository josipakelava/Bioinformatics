#include <iostream>
#include <ctime>
#include <chrono>
#include "../lib/libbf/bf/all.hpp"
#include "fasta_io.hpp"
#include "kmer_util.hpp"
#include "algorithms.hpp"

using namespace std;

template<class BF>
void doQueries(BF& bf, vector<bool>& result, const unordered_set<kmer_t>& querySet, const string& name) {
    result.resize(querySet.size());
    auto qStart = chrono::system_clock::now();
    int i = 0;
    for (kmer_t kmer : querySet) {
        result[i++] = bf.lookup(kmer);
    }
    auto qEnd = chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = qEnd - qStart;
    cout << name << " query time: " << elapsed_seconds.count() << endl;
}

void testFPRBasic(const unordered_set<kmer_t>& allKmers, const unordered_set<kmer_t>& querySet, vector<bool>& result) {
    auto basicStart = chrono::system_clock::now();
    BasicBF basicBF(allKmers);
    auto basicEnd = chrono::system_clock::now();

    std::chrono::duration<double> elapsed_seconds = basicEnd-basicStart;
    cout << "Basic init time: " << elapsed_seconds.count() << endl;
    cout << "Basic no. of kmers: " << allKmers.size() << endl;
    doQueries<BasicBF>(basicBF, result, querySet, "Basic");
}

void testOneSided(const unordered_set<kmer_t>& allKmers, const unordered_set<kmer_t>& querySet, vector<bool>& result) {
    auto start = chrono::system_clock::now();
    OneSided oneSidedBF(allKmers);
    auto end = chrono::system_clock::now();


    std::chrono::duration<double> elapsed_seconds = end-start;
    cout << "One sided init time: " << elapsed_seconds.count() << endl;

    doQueries<OneSided>(oneSidedBF, result, querySet, "1 Sided");
}

void testTwoSided(const vector<string>& sequences, const unordered_set<kmer_t>& querySet, vector<bool>& result) {
    auto start = chrono::system_clock::now();
    unordered_set<kmer_t> allKmersEdges;
    int numberOfPotentialEdges = 0;
    auto allKmers = generate_set_with_edges(sequences, allKmersEdges, numberOfPotentialEdges);
    TwoSided twoSidedBF(allKmers, allKmersEdges);
    auto end = chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end-start;
    cout << "Two sided init time: " << elapsed_seconds.count() << endl;
    cout << "Two sided potential edges: " << numberOfPotentialEdges << endl;
    cout << "Two sided edges: " << allKmersEdges.size() << endl;
    auto twoSidedEnd = chrono::system_clock::now();

    doQueries<TwoSided>(twoSidedBF, result, querySet, "2 Sided");
}

void testFPRBestFit(const vector<string>& sequences, const unordered_set<kmer_t>& querySet, vector<bool>& result) {
    unordered_set<kmer_t> bestFitKmersEdges;
    auto start = chrono::system_clock::now();
    unordered_set<kmer_t> bestFitKmers(generate_best_fit_set(sequences, bestFitKmersEdges, 1));
    SparseKBF bestFitBF(bestFitKmers, bestFitKmersEdges, 1);
    auto end = chrono::system_clock::now();

    std::chrono::duration<double> elapsed_seconds = end-start;
    cout << "Best match init time: " << elapsed_seconds.count() << endl;
    cout << "Best match no. of kmers: " << bestFitKmers.size() << endl;
    doQueries<SparseKBF>(bestFitBF, result, querySet, "Best fit");
}

void testFPRHittingSet(const vector<string>& sequences, const unordered_set<kmer_t>& querySet, vector<bool>& result) {
    unordered_set<kmer_t> hittingSetEdges;
    auto start = chrono::system_clock::now();
    unordered_set<kmer_t> hittingSetKmers(generate_hitting_set_kmers(sequences, hittingSetEdges));
    SparseRelaxedKBF hittingSetBF(hittingSetKmers, hittingSetEdges, 1);
    auto end = chrono::system_clock::now();

    std::chrono::duration<double> elapsed_seconds = end-start;
    cout << "Hit set init time: " << elapsed_seconds.count() << endl;

    cout << "Hit set no. of kmers: " << hittingSetKmers.size() << endl;
    doQueries<SparseRelaxedKBF>(hittingSetBF, result, querySet, "Hit fit");
}

void countFPR(const vector<bool>& real, const vector<bool>& result, const string& name) {
    int sum = 0;
    for(int i = 0; i < real.size(); i++) {
        sum += real[i] != result[i];
        if(real[i] && !result[i]) cout << "FALSE NEGATIVE!" << endl;
    }

    cout << name << " FPR: " << ((double) sum / real.size()) << endl;
}

void testFPR(const vector<string>& sequences) {
    unordered_set<kmer_t> allKmers(generate_set(sequences));
    unordered_set<kmer_t> querySet = query_set(allKmers);

    auto sz = querySet.size();
    vector<bool> real(sz);
    int i = 0;
    for (kmer_t kmer : querySet) {
        real[i++] = allKmers.count(kmer) == 1;
    }

    vector<bool> result;
    testFPRBasic(allKmers, querySet, result);
    countFPR(real, result, "Basic");
    testOneSided(allKmers, querySet, result);
    countFPR(real, result, "One sided");
    testTwoSided(sequences, querySet, result);
    countFPR(real, result, "Two sided");
    testFPRBestFit(sequences, querySet, result);
    countFPR(real, result, "Best fit");
}


int main(int argc, char **argv) {

    vector<string> sequences;
    read_fasta(argv[1], sequences);

    testFPR(sequences);

    return 0;
}
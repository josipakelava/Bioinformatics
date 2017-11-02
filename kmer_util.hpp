#ifndef BIOINFORMATIKA_PROJEKT_KMER_UTIL_HPP
#define BIOINFORMATIKA_PROJEKT_KMER_UTIL_HPP

#include <stdint.h>
#include <string>
#include "DEBUG_UTIL.hpp"

using namespace std;

const int KMER_LENGTH = 20;
using kmer_t = uint64_t;

const kmer_t INVALID_KMER = -1;

inline kmer_t base_to_bits(char c) {
    if (c == 'A' || c == 'a') return 0;
    if (c == 'T' || c == 't') return 1;
    if (c == 'G' || c == 'g') return 2;
    if (c == 'C' || c == 'c') return 3;

    return 0;
}

inline char bits_to_base(kmer_t kmer) {
    const char m[4] = {'A', 'T', 'G', 'C'};
    return m[kmer&3];
}

kmer_t substring_to_kmer(const string& s, int pos) {
    if(pos + KMER_LENGTH >= s.length()) {
        return INVALID_KMER;
    }

    kmer_t kmer = 0;
    for(int i = 0; i < KMER_LENGTH; i++) {
        kmer = kmer << 2 | base_to_bits(s[i]);
    }

    return kmer;
}

kmer_t string_to_kmer(const string& s) {
    return substring_to_kmer(s, 0);
}

string kmer_to_string(kmer_t kmer) {
    string s = "";
    for(int i = 0; i < KMER_LENGTH; i++) {
        s = s + bits_to_base(kmer);
        kmer >>= 2;
    }
    reverse(s.begin(), s.end());
    return s;
}

#endif //BIOINFORMATIKA_PROJEKT_KMER_UTIL_HPP

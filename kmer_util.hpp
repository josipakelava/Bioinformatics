#ifndef BIOINFORMATIKA_PROJEKT_KMER_UTIL_HPP
#define BIOINFORMATIKA_PROJEKT_KMER_UTIL_HPP

#include <stdint.h>
#include <string>
#include "DEBUG_UTIL.hpp"

using namespace std;

const int KMER_LENGTH = 20;
const int KMER_SHIFT_LEFT = (2 * (KMER_LENGTH - 1));
using kmer_t = uint64_t;

const kmer_t INVALID_KMER = -1;
const kmer_t KMER_MASK = ((uint64_t)1 << (2*KMER_LENGTH)) - 1;

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

inline kmer_t add_base_right(kmer_t kmer, char base) {
    return ((kmer << 2) & KMER_MASK) | base_to_bits(base);
}

inline kmer_t add_base_left(kmer_t kmer, char base) {
    return (kmer >> 2) | base_to_bits(base) << KMER_SHIFT_LEFT;
}

#endif //BIOINFORMATIKA_PROJEKT_KMER_UTIL_HPP

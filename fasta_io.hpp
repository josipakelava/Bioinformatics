#ifndef BIOINFORMATIKA_PROJEKT_FASTA_INPUT_H
#define BIOINFORMATIKA_PROJEKT_FASTA_INPUT_H

#include <fstream>
#include <iostream>

using namespace std;

// trim from start (in place)
static inline void ltrim(std::string &s) {
    s.erase(s.begin(), std::find_if(s.begin(), s.end(), [](int ch) {
        return !std::isspace(ch);
    }));
}

// trim from end (in place)
static inline void rtrim(std::string &s) {
    s.erase(std::find_if(s.rbegin(), s.rend(), [](int ch) {
        return !std::isspace(ch);
    }).base(), s.end());
}

// trim from both ends (in place)
static inline void trim(std::string &s) {
    ltrim(s);
    rtrim(s);
}

// trim from start (copying)
static inline std::string ltrim_copy(std::string s) {
    ltrim(s);
    return s;
}

// trim from end (copying)
static inline std::string rtrim_copy(std::string s) {
    rtrim(s);
    return s;
}

// trim from both ends (copying)
static inline std::string trim_copy(std::string s) {
    trim(s);
    return s;
}

void read_fasta(const char *filePath, vector<string>& reads) {
    std::cout << filePath << std::endl;
    ifstream infile(filePath);
    string buildSeq = "";
    for(string line; getline(infile, line);) {
        line = trim_copy(line);
        if(line.empty()) continue;

        if(line[0] == '>') {
            if(!buildSeq.empty()) {
                reads.emplace_back(buildSeq);
            }
            buildSeq = "";
        } else {
            buildSeq += line;
        }
    }

    if(!buildSeq.empty()) reads.emplace_back(buildSeq);
}

#endif //BIOINFORMATIKA_PROJEKT_FASTA_INPUT_H
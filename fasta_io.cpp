//
// Created by Ivan Krpelnik on 02/11/17.
//

#include "fasta_io.h"

#include <fstream>
#include <iostream>
#include <boost/algorithm/string.hpp>

void read_fasta(const char *filePath, vector<string>& reads) {
    std::cout << filePath << std::endl;
    ifstream infile(filePath);
    string buildSeq = "";
    for(string line; getline(infile, line);) {
        line = boost::trim_copy(line);
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
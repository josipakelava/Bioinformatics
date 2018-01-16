cmake_minimum_required(VERSION 3.8)
project(Bioinformatika_projekt)

set(CMAKE_CXX_STANDARD 17)

set(SOURCE_FILES
        Project/main.cpp
        lib/libbf/src/bitvector.cpp
        lib/libbf/src/counter_vector.cpp
        lib/libbf/src/hash.cpp
        lib/libbf/src/bloom_filter/a2.cpp
        lib/libbf/src/bloom_filter/basic.cpp
        lib/libbf/src/bloom_filter/bitwise.cpp
        lib/libbf/src/bloom_filter/counting.cpp
        lib/libbf/src/bloom_filter/stable.cpp
        Project/kmer_util.hpp
        Project/DEBUG_UTIL.hpp
        Project/algorithms.hpp
        Project/fasta_io.hpp)

add_executable(Bioinformatika_projekt ${SOURCE_FILES})

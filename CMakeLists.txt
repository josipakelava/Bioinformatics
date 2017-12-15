cmake_minimum_required(VERSION 3.8)
project(Bioinformatika_projekt)

set(CMAKE_CXX_STANDARD 17)

set(SOURCE_FILES
        main.cpp
        lib/libbf/src/bitvector.cpp
        lib/libbf/src/counter_vector.cpp
        lib/libbf/src/hash.cpp
        lib/libbf/src/bloom_filter/a2.cpp
        lib/libbf/src/bloom_filter/basic.cpp
        lib/libbf/src/bloom_filter/bitwise.cpp
        lib/libbf/src/bloom_filter/counting.cpp
        lib/libbf/src/bloom_filter/stable.cpp
        kmer_util.hpp
        DEBUG_UTIL.hpp
        algorithms.hpp
        fasta_io.hpp)

find_package(Boost COMPONENTS system filesystem REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})

add_executable(Bioinformatika_projekt ${SOURCE_FILES})
target_link_libraries(Bioinformatika_projekt ${Boost_LIBRARIES})
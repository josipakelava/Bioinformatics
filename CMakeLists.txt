cmake_minimum_required(VERSION 3.8)
project(Bioinformatika_projekt)

set(CMAKE_CXX_STANDARD 17)

find_package(Boost COMPONENTS system filesystem REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})

find_package(libbf)

set(SOURCE_FILES main.cpp fasta_io.cpp fasta_io.h kmer_util.hpp DEBUG_UTIL.hpp)
add_executable(Bioinformatika_projekt ${SOURCE_FILES})
target_link_libraries(Bioinformatika_projekt ${Boost_LIBRARIES})
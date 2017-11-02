cmake_minimum_required(VERSION 3.8)
project(Bioinformatika_projekt)

set(CMAKE_CXX_STANDARD 17)

find_package(Boost COMPONENTS system filesystem REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})

set(SOURCE_FILES main.cpp fasta_io.cpp fasta_io.h)
add_executable(Bioinformatika_projekt ${SOURCE_FILES})
target_link_libraries(Bioinformatika_projekt ${Boost_LIBRARIES})
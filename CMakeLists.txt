cmake_minimum_required(VERSION 3.5)
project(SWPraktikum)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")


# Search for zlib as a dependency for SeqAn.
find_package (ZLIB)

# Load the SeqAn module and fail if not found.
find_package (SeqAn REQUIRED)

# Add include directories.
include_directories (${SEQAN_INCLUDE_DIRS})

# Add definitions set by find_package (SeqAn).
add_definitions (${SEQAN_DEFINITIONS})







set(SOURCE_FILES main.cpp)
add_executable(SWPraktikum ${SOURCE_FILES})

target_link_libraries (SWPraktikum ${SEQAN_LIBRARIES})
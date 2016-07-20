cmake_minimum_required (VERSION 2.8)
project (hts_common)

file(GLOB SOURCES "src/*.cpp")

add_library(hts_common SHARED ${SOURCES})


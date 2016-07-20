cmake_minimum_required (VERSION 2.8)
project (hts_common)

file(GLOB SOURCES "*.cpp")

add_library(hts_common SHARED ${SOURCES})


cmake_minimum_required (VERSION 2.8)
project (HTStream)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

if (${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")
endif()
if (${CMAKE_CXX_COMPILER_ID} STREQUAL "Clang")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Weverything -Wno-c++98-compat-pedantic -Wno-padded")
endif()

message( status "cmake cxx flags: " ${CMAKE_CXX_FLAGS})

include(gtest.cmake)

enable_testing()

add_subdirectory(common)
add_subdirectory(Super-Deduper)

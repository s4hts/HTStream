cmake_minimum_required (VERSION 2.8)
project (HTStream)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

if (${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")
endif()
if (${CMAKE_CXX_COMPILER_ID} STREQUAL "Clang")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Weverything -Wno-c++98-compat-pedantic -Wno-padded -Wno-undef -Wno-weak-vtables -Wno-global-constructors -Wno-missing-prototypes")
## todo libc++ 
# -stdlib=libc++")
#  set(PROJECT_LINK_LIBS "-lc++ -lc++abi -lm -lc -lgcc_s -lgcc")
endif()

message( status "cmake cxx flags: " ${CMAKE_CXX_FLAGS})

include(gtest.cmake)
message( status "include dirs: " ${GTEST_INCLUDE_DIRS})

enable_testing()

add_subdirectory(common)
add_subdirectory(Super-Deduper)
add_subdirectory(Tab-Converter)

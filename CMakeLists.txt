cmake_minimum_required (VERSION 3.2)
project (HTStream)

option(BUILD_STATIC_BIN "builds static binary (linux only)" OFF)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 ")

set(Boost_USE_STATIC_LIBS OFF)
set(BUILD_SHARED_LIBS ON)

if (APPLE)
elseif (UNIX)
  if (BUILD_STATIC_BIN)
    set(EXTRA_BUILD_FLAGS "-static")
    set(Boost_USE_STATIC_LIBS ON)
    set(BUILD_SHARED_LIBS OFF)
  endif()
endif()

if (${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")
endif()
if (${CMAKE_CXX_COMPILER_ID} STREQUAL "Clang")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Weverything -Wno-c++98-compat-pedantic -Wno-padded -Wno-undef -Wno-weak-vtables -Wno-global-constructors -Wno-missing-prototypes -Wundefined-func-template")
## todo libc++
# -stdlib=libc++")
#  set(PROJECT_LINK_LIBS "-lc++ -lc++abi -lm -lc -lgcc_s -lgcc")
endif()

message( status "cmake cxx flags: " ${CMAKE_CXX_FLAGS})

include(gtest.cmake)
message( status "include dirs: " ${GTEST_INCLUDE_DIRS})

enable_testing()

add_subdirectory(common)
add_subdirectory(hts_Stats)
add_subdirectory(hts_SuperDeduper)
add_subdirectory(hts_PolyATTrim)
add_subdirectory(hts_QWindowTrim)
add_subdirectory(hts_CutTrim)
add_subdirectory(hts_SeqScreener)
add_subdirectory(hts_NTrimmer)
add_subdirectory(hts_AdapterTrimmer)
add_subdirectory(hts_Overlapper)

cmake_minimum_required (VERSION 2.8)
project (Super-Deduper)

file(GLOB SOURCES "source/*.cpp")

include_directories(source)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

if (${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")
endif()
if (${CMAKE_CXX_COMPILER_ID} STREQUAL "Clang")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Weverything -Wno-c++98-compat-pedantic -Wno-padded")
endif()

message( status "cmake cxx flags: " ${CMAKE_CXX_FLAGS})

add_executable(super_deduper ${SOURCES})

add_custom_target(test
  DEPENDS super_deduper
  COMMAND python ${CMAKE_SOURCE_DIR}/test_suite/testSD.py # todo add argument for super d
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/test_suite
)

add_custom_target(valgrind
  DEPENDS super_deduper
  COMMAND valgrind ./super_deduper -F -1 ${CMAKE_SOURCE_DIR}/test_suite/fastqFiles/testCase_1X_R1.fastq -2 ${CMAKE_SOURCE_DIR}/test_suite/fastqFiles/testCase_1X_R1.fastq -O
)

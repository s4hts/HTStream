cmake_minimum_required (VERSION 3.2)
project (Super-Deduper)

file(GLOB SOURCES "*.cpp")

include_directories(${CMAKE_SOURCE_DIR}/common)

set(PROJECT_LINK_LIBS libhts_common.so)
link_directories(${CMAKE_BINARY_DIR}/common)

add_executable(super_deduper ${SOURCES})
target_link_libraries(super_deduper ${PROJECT_LINK_LIBS})

add_custom_target(testSD
  DEPENDS super_deduper
  COMMAND ${CMAKE_COMMAND} -E env "PATH=$ENV{PATH}:${CMAKE_BINARY_DIR}/Super-Deduper" python ${CMAKE_SOURCE_DIR}/Super-Deduper/test_suite/testSD.py
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/Super-Deduper/test_suite
)

add_custom_target(valgrind
  DEPENDS super_deduper
  COMMAND valgrind ./super_deduper -F -1 ${CMAKE_SOURCE_DIR}/test_suite/fastqFiles/testCase_1X_R1.fastq -2 ${CMAKE_SOURCE_DIR}/test_suite/fastqFiles/testCase_1X_R1.fastq -O
)

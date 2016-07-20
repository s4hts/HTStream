cmake_minimum_required (VERSION 3.2)
project (Super-Deduper)

file(GLOB SOURCES "src/*.cpp")

set(COMMON_INCLUDES ${CMAKE_SOURCE_DIR}/common/src)
include_directories(${COMMON_INCLUDES})

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

## test
set(PROJECT_TEST_NAME ${PROJECT_NAME}_test)

include_directories(${GTEST_INCLUDE_DIRS} ${COMMON_INCLUDES})

file(GLOB TEST_SRC_FILES ${PROJECT_SOURCE_DIR}/test/*.cpp)
message( status "inlcude dirs: " ${GTEST_INCLUDE_DIRS})

add_executable(${PROJECT_TEST_NAME} ${TEST_SRC_FILES})
add_dependencies(${PROJECT_TEST_NAME} googletest)
add_test(test_superd ${PROJECT_TEST_NAME})

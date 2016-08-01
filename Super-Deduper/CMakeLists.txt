cmake_minimum_required (VERSION 3.2)
project (Super-Deduper)

file(GLOB sources "src/*.cpp")

message( status "sources: " ${sources})


set(COMMON_INCLUDES ${CMAKE_SOURCE_DIR}/common/src)

set(PROJECT_LINK_LIBS hts_common)
link_directories(${CMAKE_BINARY_DIR}/common)

FIND_PACKAGE( Boost 1.54 COMPONENTS program_options REQUIRED )
include_directories(${COMMON_INCLUDES} ${Boost_INCLUDE_DIR} )


add_executable(super_deduper ${sources})
target_link_libraries(super_deduper ${PROJECT_LINK_LIBS} ${Boost_LIBRARIES})

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
find_package(Threads REQUIRED)
include_directories(${GTEST_INCLUDE_DIRS} ${COMMON_INCLUDES} ${CMAKE_SOURCE_DIR}/Super-Deduper/src)

file(GLOB TEST_SRC_FILES ${PROJECT_SOURCE_DIR}/test/*.cpp)
link_directories(${GTEST_LIBS_DIR})
list(REMOVE_ITEM sources ${CMAKE_CURRENT_SOURCE_DIR}/src/super-deduper.cpp)
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Dprivate=public")
add_executable(${PROJECT_TEST_NAME} ${TEST_SRC_FILES} ${sources})
#add_dependencies(${PROJECT_TEST_NAME} googletest)
target_link_libraries(${PROJECT_TEST_NAME} gtest gtest_main ${CMAKE_THREAD_LIBS_INIT} ${PROJECT_LINK_LIBS})
add_test(test_superd ${PROJECT_TEST_NAME})

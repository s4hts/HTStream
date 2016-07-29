cmake_minimum_required (VERSION 3.2)
project (hts_common)

file(GLOB sources "src/*.cpp")
message( status "sources: " ${sources})

set(COMMON_INCLUDES ${CMAKE_SOURCE_DIR}/common/src)
add_library(hts_common SHARED ${sources})


## test
set(PROJECT_TEST_NAME ${PROJECT_NAME}_test)
find_package(Threads REQUIRED)
include_directories(${GTEST_INCLUDE_DIRS} ${COMMON_INCLUDES} )

file(GLOB TEST_SRC_FILES ${PROJECT_SOURCE_DIR}/test/*.cpp)
link_directories(${GTEST_LIBS_DIR})

add_executable(${PROJECT_TEST_NAME} ${TEST_SRC_FILES} ${sources})
#add_dependencies(${PROJECT_TEST_NAME} googletest)
target_link_libraries(${PROJECT_TEST_NAME} gtest gtest_main ${CMAKE_THREAD_LIBS_INIT} ${PROJECT_LINK_LIBS})
add_test(test_common ${PROJECT_TEST_NAME})



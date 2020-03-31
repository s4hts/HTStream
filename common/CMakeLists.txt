cmake_minimum_required (VERSION 3.2)
project (hts_common)

add_custom_command(
  OUTPUT 
    ${CMAKE_BINARY_DIR}/common/version.h
  COMMAND ${CMAKE_COMMAND} -D SRC_DIR=${CMAKE_SOURCE_DIR} -D BIN_DIR=${CMAKE_BINARY_DIR} -P ${CMAKE_SOURCE_DIR}/cmake/version.cmake
  DEPENDS always_rebuild
  COMMENT "Configuring version.h")

add_custom_target(generateVersion ALL
  DEPENDS ${CMAKE_BINARY_DIR}/common/version.h
  COMMENT "generateVersion target")

# Dummy output which is never actually produced. Anything that depends on
# this will always be rebuilt.
add_custom_command(
    OUTPUT always_rebuild
    COMMAND echo -n
    )  

file(GLOB sources "src/*.cpp")
list(APPEND sources  ${CMAKE_BINARY_DIR}/common/version.h)
message( "common sources: " ${sources})

FIND_PACKAGE( Boost 1.56 COMPONENTS system program_options filesystem iostreams REQUIRED )

include_directories(${GTEST_INCLUDE_DIRS} ${Boost_INCLUDE_DIR} ${COMMON_INCLUDES} )

add_library(${PROJECT_NAME} STATIC ${sources})
target_link_libraries(${PROJECT_NAME}  ${Boost_LIBRARIES})

## test
set(PROJECT_TEST_NAME ${PROJECT_NAME}_test)
find_package(Threads REQUIRED)

include_directories(${GTEST_INCLUDE_DIRS} ${Boost_INCLUDE_DIR} ${COMMON_INCLUDES} )
file(GLOB TEST_SRC_FILES ${PROJECT_SOURCE_DIR}/test/*.cpp)
link_directories(${GTEST_LIBS_DIR})

add_executable(${PROJECT_TEST_NAME} ${TEST_SRC_FILES} ${sources})

#add_dependencies(${PROJECT_TEST_NAME} googletest)
target_link_libraries(${PROJECT_TEST_NAME} gtest gtest_main ${CMAKE_THREAD_LIBS_INIT} ${Boost_LIBRARIES} ${PROJECT_LINK_LIBS})
add_test(test_common ${PROJECT_TEST_NAME})

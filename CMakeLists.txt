cmake_minimum_required (VERSION 3.2)
project (HTStream)

option(BUILD_STATIC_BIN "builds static binary (linux only)" OFF)
option(BUILD_GPERFTOOLS "links to google perftools" OFF)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake")
message( status "cmake module path: " ${CMAKE_MODULE_PATH})

set(CMAKE_CXX_STANDARD 11)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(Boost_USE_STATIC_LIBS OFF)
set(BUILD_SHARED_LIBS ON)

if (APPLE)
elseif (UNIX)
  if (BUILD_STATIC_BIN)
    set(EXTRA_BUILD_FLAGS ${EXTRA_BUILD_FLAGS} "-static -pthread -Wl,--whole-archive -lpthread -Wl,--no-whole-archive")
    set(Boost_USE_STATIC_LIBS ON)
    set(BUILD_SHARED_LIBS OFF)
  endif()
endif()

if (${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Werror")
endif()
if (${CMAKE_CXX_COMPILER_ID} STREQUAL "Clang")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Werror -Wextra")
## todo libc++
# -stdlib=libc++")
#  set(PROJECT_LINK_LIBS "-lc++ -lc++abi -lm -lc -lgcc_s -lgcc")
endif()

message( status "cmake cxx flags: " ${CMAKE_CXX_FLAGS})

include(cmake/gtest.cmake)
message( status "include dirs: " ${GTEST_INCLUDE_DIRS})

enable_testing()

set(COMMON_INCLUDES ${CMAKE_SOURCE_DIR}/common/src ${CMAKE_BINARY_DIR}/common)

add_subdirectory(common)
add_subdirectory(hts_AdapterTrimmer)
add_subdirectory(hts_CutTrim)
add_subdirectory(hts_LengthFilter)
add_subdirectory(hts_NTrimmer)
add_subdirectory(hts_Overlapper)
add_subdirectory(hts_PolyATTrim)
add_subdirectory(hts_Primers)
add_subdirectory(hts_QWindowTrim)
add_subdirectory(hts_SeqScreener)
add_subdirectory(hts_Stats)
add_subdirectory(hts_SuperDeduper)

add_custom_target( copy-compile-commands ALL
    COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CMAKE_BINARY_DIR}/compile_commands.json ${CMAKE_CURRENT_LIST_DIR}
    COMMENT "copy compile_commands"
    )

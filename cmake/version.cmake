cmake_minimum_required (VERSION 3.2)

execute_process(
  COMMAND git describe --tags --dirty
  OUTPUT_VARIABLE HTSTREAM_VERSION
  OUTPUT_STRIP_TRAILING_WHITESPACE
  RESULT_VARIABLE result
)

if(result EQUAL "0")
  message("execute git version " ${HTSTREAM_VERSION})
  configure_file(${SRC_DIR}/common/src/version.h.in
    ${BIN_DIR}/common/version.h @ONLY)
else()
  message("using release version")
  configure_file(${SRC_DIR}/common/src/version.h.release
    ${BIN_DIR}/common/version.h COPYONLY)
endif()

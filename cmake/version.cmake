cmake_minimum_required (VERSION 3.2)

execute_process(
  COMMAND git describe --tags
  OUTPUT_VARIABLE HTSTREAM_VERSION
  OUTPUT_STRIP_TRAILING_WHITESPACE
)
message("execute git version " ${HTSTREAM_VERSION})


configure_file(${SRC_DIR}/common/src/version.h.in
  ${BIN_DIR}/common/version.h @ONLY)

# Repository Guidelines

## Project Structure & Module Organization
- Core shared code in `common/src` with unit tests in `common/test`; executables live in tool-specific folders (`hts_AdapterTrimmer`, `hts_SuperDeduper`, etc.), each with `src/`, `test/`, and a local `README` describing options.
- CMake configuration lives in the repo root and `cmake/`; third-party archives (googletest) are in `ext/`.
- `regression/` holds small FASTQ fixtures and expected outputs for end-to-end checks. Build artifacts should stay outside the repo (e.g., `build/`).

## Build, Test, and Development Commands
- Configure an out-of-source build: `cmake -S . -B build -DCMAKE_BUILD_TYPE=Release` (add `-DBUILD_STATIC_BIN=ON` for static Linux binaries or `-DBUILD_GPERFTOOLS=ON` when profiling).
- Compile: `cmake --build build -j$(nproc)`; copy of `compile_commands.json` is generated for tooling.
- Unit tests: `ctest --test-dir build` (or `cmake --build build --target test` / `make test` inside `build`).
- Regression check after building: `cd regression && ./regression.sh ../build test` to compare outputs against fixtures.

## Coding Style & Naming Conventions
- C++11 codebase compiled with `-Wall -Werror` (and `-Wextra` on Clang); keep code warning-clean.
- Prefer 4-space indentation, Allman-style braces, and explicit `std::`/`boost::` qualifiers; headers guard with `#ifndef/#define`.
- Executable targets follow `hts_*` naming; tests mirror module names (e.g., `hts_TestQWindowTrim.cpp`).
- Use existing option-parsing helpers in `common/src/utils.h` and keep new flags documented in the module `README`.

## Testing Guidelines
- Unit tests use GoogleTest (vendored under `ext/googletest`); place new cases in the module’s `test/` folder and keep filenames consistent with the tool under test.
- When behavior depends on streaming order, add deterministic sorting (see `regression/regression.sh` handling for `hts_SuperDeduper`).
- For new CLI flags, include tests covering default behavior, error handling (invalid parameters), and representative FASTQ inputs.

## Commit & Pull Request Guidelines
- Commit messages are short and imperative (e.g., “Fix clang warnings”, “Update release doc”); reference issues with `#123` when relevant.
- Before opening a PR, run unit tests and regression checks, and update version info (`common/src/version.h.release`) only when preparing a release branch.
- PRs should summarize the change, outline user-visible flags or defaults that shifted, note test coverage, and link any tracking issue. Include sample commands if the workflow changes.

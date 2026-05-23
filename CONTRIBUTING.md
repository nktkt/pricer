# Contributing to pricer

Thanks for your interest in **pricer** — a header-only C++17 quant library with an
LLVM-JIT payoff compiler. This guide covers how to build, test, and add to the
project. The core stays small, heavily commented, and dependency-light; please
keep new work in that spirit.

## Prerequisites

- **CMake 3.16+** and a **C++17 compiler** (`clang++` or `g++`).
- **LLVM** *(optional)* — only needed for the JIT examples (`jit_payoff`,
  `path_dependent`, `simd_payoff`) and the `test_payoff_dsl` test. On macOS:
  ```sh
  brew install llvm
  ```
  If CMake doesn't find it automatically, point it there:
  ```sh
  cmake -S . -B build -DLLVM_DIR=$(llvm-config --cmakedir)
  ```
- **Python 3.9+** *(optional)* — for the pybind11 bindings (`pip install .`).
- **Docker** *(optional)* — for building and running the REST server container.

## Build & test

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel
ctest --test-dir build --output-on-failure
```

The JIT components are **ON by default** but are only built when LLVM is found.
CI builds with `-DPRICER_ENABLE_JIT=OFF` (the hosted runners ship an LLVM whose
ORC JIT cannot initialize there), so the JIT path is exercised **locally**. There
are **30 tests locally** (JIT on) and **29 in CI** (JIT off).

### CMake options

| Option | Default | Purpose |
|--------|---------|---------|
| `PRICER_BUILD_EXAMPLES` | `ON` | Build the example programs |
| `PRICER_BUILD_TESTS` | `ON` | Build the CTest suite |
| `PRICER_ENABLE_JIT` | `ON` | Build the LLVM JIT components (needs a working LLVM) |
| `PRICER_BUILD_PYTHON` | `OFF` | Build the Python bindings (needs pybind11) |
| `PRICER_BUILD_SERVER` | `OFF` | Build the REST pricing service (POSIX sockets) |

## Project conventions

- **The core is header-only.** The `pricer` library is an `INTERFACE` CMake
  target — add new functionality as headers under `include/pricer/`.
- **First-party targets are warning-clean.** Examples, tests, and the server link
  the `pricer_warnings` interface target, which compiles with `-Wall -Wextra
  -Werror`. Keep the build free of warnings.
- **Tests are dependency-free.** Each test is a standalone executable that uses
  `tests/check.hpp`, a tiny assert helper (`check::approx`, `check::is_true`,
  `check::report`). Each returns a non-zero exit code on failure and is registered
  with CTest.
- **Pair features with coverage.** Every feature should come with a test (in
  `tests/`) and, ideally, a runnable example (in `examples/`).

## Adding a feature

1. Add the header to `include/pricer/`.
2. Add a `test_<name>.cpp` to `tests/`, then append `test_<name>` to the
   `foreach` list in `tests/CMakeLists.txt`.
3. *(Optional but encouraged)* add a `<name>_demo.cpp` to `examples/` and append
   it to the `foreach` list in `examples/CMakeLists.txt`.
4. Add a `### Added` bullet under `## [Unreleased]` in `CHANGELOG.md`.
5. Build with JIT both **on** and **off**, and run `ctest`.

## Code style

- **C++17.** Match the surrounding code's naming, comment density, and idiom.
- Each header begins with a short comment explaining **what** it does and **why**.
- Prefer the existing patterns over introducing new dependencies — the library
  core has zero heavy dependencies, and keeping it that way is a feature.

## Commits & pull requests

The project follows [Keep a Changelog](https://keepachangelog.com/) and
[Semantic Versioning](https://semver.org/) (see [`CHANGELOG.md`](CHANGELOG.md)).

- Write focused commits with clear messages.
- Update the `[Unreleased]` section of `CHANGELOG.md`.
- Ensure CI is green.

CI runs on **Linux and macOS** and consists of:

- a **C++ build/test matrix** (`build-and-test`),
- a **`python`** job (`pip install .` + the Python tests),
- a **`server`** job (builds the REST server + runs `server/smoke_test.py`),
- a **`docker`** job (builds the image + checks the container endpoints),
- plus a separate **`docs`** workflow (Doxygen → GitHub Pages).

## Release process (maintainers)

1. Bump the version in `CMakeLists.txt`, `Doxyfile`, `pyproject.toml`, and
   `python/src/bindings.cpp`.
2. In `CHANGELOG.md`, promote `## [Unreleased]` to a new `## [X.Y.Z]` section with
   the date, and update the link references at the bottom.
3. Verify CI is green on the release commit.
4. Tag and publish:
   ```sh
   git tag -a vX.Y.Z
   git push origin vX.Y.Z
   gh release create vX.Y.Z --verify-tag --notes-file <notes>
   ```

## See also

- [`README.md`](README.md) — usage and examples
- [`ARCHITECTURE.md`](ARCHITECTURE.md) — design
- [`ROADMAP.md`](ROADMAP.md) — the plan
- **API reference:** https://nktkt.github.io/pricer/

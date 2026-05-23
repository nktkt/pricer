cmake_minimum_required(VERSION 3.16)
# C is enabled because LLVM's CMake package runs C-based feature checks.
project(pricer LANGUAGES C CXX VERSION 0.12.0)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "" FORCE)
endif()

# Header-only core library.
add_library(pricer INTERFACE)
add_library(pricer::pricer ALIAS pricer)
target_include_directories(pricer INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>)
target_compile_features(pricer INTERFACE cxx_std_17)

# Strict warnings, treated as errors, for first-party targets. Linked into our
# own examples/tests/server (not the LLVM- or pybind-heavy targets, whose
# third-party headers we do not control).
add_library(pricer_warnings INTERFACE)
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
  target_compile_options(pricer_warnings INTERFACE -Wall -Wextra -Werror)
endif()
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  # The 256-bit SIMD vector types (simd.hpp) trip GCC's -Wpsabi note about the
  # GCC-4.6 ABI change for >16-byte vectors passed by value. Those types are
  # always inlined and never cross a real ABI boundary here, so the note is
  # spurious; silence it (it would otherwise be fatal under -Werror).
  target_compile_options(pricer_warnings INTERFACE -Wno-psabi)
endif()

option(PRICER_BUILD_EXAMPLES "Build example programs" ON)
option(PRICER_BUILD_TESTS "Build tests" ON)
# LLVM JIT components (payoff DSL). ON by default, but only built when LLVM is
# also found. Turn OFF where a discoverable LLVM cannot actually host a JIT
# (e.g. some CI runners) so the rest of the project still builds and tests.
option(PRICER_ENABLE_JIT "Build LLVM JIT components (needs a working LLVM)" ON)
option(PRICER_BUILD_PYTHON "Build the Python bindings (needs pybind11)" OFF)
option(PRICER_BUILD_SERVER "Build the REST pricing service (POSIX sockets; Linux/macOS)" OFF)

if(PRICER_BUILD_EXAMPLES)
  add_subdirectory(examples)
endif()

if(PRICER_BUILD_TESTS)
  enable_testing()
  add_subdirectory(tests)
endif()

if(PRICER_BUILD_PYTHON)
  add_subdirectory(python)
endif()

if(PRICER_BUILD_SERVER)
  add_subdirectory(server)
endif()

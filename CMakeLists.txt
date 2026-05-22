cmake_minimum_required(VERSION 3.16)
# C is enabled because LLVM's CMake package runs C-based feature checks.
project(pricer LANGUAGES C CXX VERSION 0.1.0)

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

option(PRICER_BUILD_EXAMPLES "Build example programs" ON)
option(PRICER_BUILD_TESTS "Build tests" ON)
# LLVM JIT components (payoff DSL). ON by default, but only built when LLVM is
# also found. Turn OFF where a discoverable LLVM cannot actually host a JIT
# (e.g. some CI runners) so the rest of the project still builds and tests.
option(PRICER_ENABLE_JIT "Build LLVM JIT components (needs a working LLVM)" ON)
option(PRICER_BUILD_PYTHON "Build the Python bindings (needs pybind11)" OFF)

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

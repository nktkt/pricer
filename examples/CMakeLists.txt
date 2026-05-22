find_package(Threads REQUIRED)

# Examples that only need the header-only core.
foreach(ex black_scholes_demo convergence parallel_mc greeks barrier_option)
  add_executable(${ex} ${ex}.cpp)
  target_link_libraries(${ex} PRIVATE pricer Threads::Threads)
endforeach()

# The LLVM JIT example is optional: built only when LLVM is available.
find_package(LLVM CONFIG QUIET)
if(LLVM_FOUND)
  message(STATUS "LLVM ${LLVM_PACKAGE_VERSION} found — building jit_payoff")
  add_executable(jit_payoff jit_payoff.cpp)
  target_include_directories(jit_payoff PRIVATE ${LLVM_INCLUDE_DIRS})
  target_compile_definitions(jit_payoff PRIVATE ${LLVM_DEFINITIONS})
  llvm_map_components_to_libnames(LLVM_LIBS core orcjit native)
  target_link_libraries(jit_payoff PRIVATE pricer ${LLVM_LIBS})
  # The DSL parser uses exceptions; LLVM may ship built with -fno-exceptions.
  target_compile_options(jit_payoff PRIVATE -fexceptions)
else()
  message(STATUS "LLVM not found — skipping jit_payoff example")
endif()

find_package(Threads REQUIRED)

# Examples that only need the header-only core.
foreach(ex black_scholes_demo convergence parallel_mc greeks barrier_option benchmark
          variance_reduction mc_greeks risk_demo pricer_cli vol_surface_demo)
  add_executable(${ex} ${ex}.cpp)
  target_link_libraries(${ex} PRIVATE pricer Threads::Threads)
endforeach()

# The LLVM JIT examples are optional: built only when JIT is enabled AND LLVM is found.
if(PRICER_ENABLE_JIT)
  find_package(LLVM CONFIG QUIET)
endif()
if(PRICER_ENABLE_JIT AND LLVM_FOUND)
  message(STATUS "LLVM ${LLVM_PACKAGE_VERSION} found — building JIT examples")
  llvm_map_components_to_libnames(LLVM_LIBS core orcjit native)
  foreach(ex jit_payoff path_dependent simd_payoff)
    add_executable(${ex} ${ex}.cpp)
    target_include_directories(${ex} PRIVATE ${LLVM_INCLUDE_DIRS})
    target_compile_definitions(${ex} PRIVATE ${LLVM_DEFINITIONS})
    target_link_libraries(${ex} PRIVATE pricer ${LLVM_LIBS})
    # The DSL parser uses exceptions; LLVM may ship built with -fno-exceptions.
    target_compile_options(${ex} PRIVATE -fexceptions)
  endforeach()
elseif(NOT PRICER_ENABLE_JIT)
  message(STATUS "PRICER_ENABLE_JIT=OFF — skipping JIT examples")
else()
  message(STATUS "LLVM not found — skipping JIT examples")
endif()

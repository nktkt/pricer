# Each test is a standalone executable registered with CTest.
find_package(Threads REQUIRED)
foreach(t test_black_scholes test_greeks test_monte_carlo test_parallel test_variance_reduction
          test_implied_vol test_greeks_mc test_risk)
  add_executable(${t} ${t}.cpp)
  target_link_libraries(${t} PRIVATE pricer Threads::Threads)
  target_include_directories(${t} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
  add_test(NAME ${t} COMMAND ${t})
endforeach()

# Smoke-test the CLI tool (only when the examples, hence the CLI, were built).
if(TARGET pricer_cli)
  add_test(NAME cli_price
           COMMAND pricer_cli price --type call --S 100 --K 100 --r 0.05 --sigma 0.2 --T 1)
  set_tests_properties(cli_price PROPERTIES PASS_REGULAR_EXPRESSION "price: 10\\.45")
endif()

# DSL/JIT test runs only when JIT is enabled and LLVM is found.
if(PRICER_ENABLE_JIT)
  find_package(LLVM CONFIG QUIET)
endif()
if(PRICER_ENABLE_JIT AND LLVM_FOUND)
  llvm_map_components_to_libnames(LLVM_LIBS core orcjit native)
  add_executable(test_payoff_dsl test_payoff_dsl.cpp)
  target_include_directories(test_payoff_dsl PRIVATE ${CMAKE_CURRENT_SOURCE_DIR} ${LLVM_INCLUDE_DIRS})
  target_compile_definitions(test_payoff_dsl PRIVATE ${LLVM_DEFINITIONS})
  target_link_libraries(test_payoff_dsl PRIVATE pricer ${LLVM_LIBS})
  target_compile_options(test_payoff_dsl PRIVATE -fexceptions)
  add_test(NAME test_payoff_dsl COMMAND test_payoff_dsl)
endif()

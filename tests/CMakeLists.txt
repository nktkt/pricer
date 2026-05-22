# Each test is a standalone executable registered with CTest.
foreach(t test_black_scholes test_greeks test_monte_carlo)
  add_executable(${t} ${t}.cpp)
  target_link_libraries(${t} PRIVATE pricer)
  target_include_directories(${t} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
  add_test(NAME ${t} COMMAND ${t})
endforeach()

# DSL/JIT test runs only when LLVM is available.
find_package(LLVM CONFIG QUIET)
if(LLVM_FOUND)
  llvm_map_components_to_libnames(LLVM_LIBS core orcjit native)
  add_executable(test_payoff_dsl test_payoff_dsl.cpp)
  target_include_directories(test_payoff_dsl PRIVATE ${CMAKE_CURRENT_SOURCE_DIR} ${LLVM_INCLUDE_DIRS})
  target_compile_definitions(test_payoff_dsl PRIVATE ${LLVM_DEFINITIONS})
  target_link_libraries(test_payoff_dsl PRIVATE pricer ${LLVM_LIBS})
  target_compile_options(test_payoff_dsl PRIVATE -fexceptions)
  add_test(NAME test_payoff_dsl COMMAND test_payoff_dsl)
endif()

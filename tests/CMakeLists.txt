# Each test is a standalone executable registered with CTest.
foreach(t test_black_scholes test_greeks test_monte_carlo)
  add_executable(${t} ${t}.cpp)
  target_link_libraries(${t} PRIVATE pricer)
  target_include_directories(${t} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
  add_test(NAME ${t} COMMAND ${t})
endforeach()

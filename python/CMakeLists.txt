# Python bindings (built when PRICER_BUILD_PYTHON=ON; needs pybind11).
find_package(pybind11 CONFIG REQUIRED)
find_package(Threads REQUIRED)

pybind11_add_module(_pricer src/bindings.cpp)
target_link_libraries(_pricer PRIVATE pricer Threads::Threads)

# For local use, drop the extension next to the package __init__.py so that
# `PYTHONPATH=python python -c "import pricer"` works after a plain CMake build.
set_target_properties(_pricer PROPERTIES
  LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/pricer)

# For wheels (scikit-build-core), install the extension into the package.
install(TARGETS _pricer LIBRARY DESTINATION pricer)

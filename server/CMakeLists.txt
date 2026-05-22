# REST pricing service (built when PRICER_BUILD_SERVER=ON). POSIX sockets only,
# so it is Linux/macOS only and has no third-party dependencies.
find_package(Threads REQUIRED)
add_executable(pricer_server pricer_server.cpp)
target_link_libraries(pricer_server PRIVATE pricer Threads::Threads)

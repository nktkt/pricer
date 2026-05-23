# Multi-stage build for the pricer REST service.
#
# Stage 1 compiles the dependency-free server; stage 2 ships only the resulting
# binary on a slim runtime, run as a non-root user. The server has no third-party
# dependencies (POSIX sockets only), so the runtime image needs just the C++
# standard library.
#
#   docker build -t pricer-server .
#   docker run --rm -p 8080:8080 pricer-server
#   curl 'http://127.0.0.1:8080/price?type=call&S=100&K=100&r=0.05&sigma=0.2&T=1'
#   curl http://127.0.0.1:8080/metrics

# --- build stage ---
FROM debian:bookworm-slim AS build
RUN apt-get update \
    && apt-get install -y --no-install-recommends build-essential cmake \
    && rm -rf /var/lib/apt/lists/*
WORKDIR /src
# Only the pieces the server target needs (examples/tests/python are off).
COPY CMakeLists.txt ./
COPY include ./include
COPY server ./server
RUN cmake -S . -B build \
        -DCMAKE_BUILD_TYPE=Release \
        -DPRICER_BUILD_SERVER=ON \
        -DPRICER_BUILD_EXAMPLES=OFF \
        -DPRICER_BUILD_TESTS=OFF \
        -DPRICER_ENABLE_JIT=OFF \
    && cmake --build build --parallel --target pricer_server

# --- runtime stage ---
FROM debian:bookworm-slim AS runtime
RUN apt-get update \
    && apt-get install -y --no-install-recommends libstdc++6 \
    && rm -rf /var/lib/apt/lists/* \
    && useradd --system --no-create-home --uid 10001 pricer
COPY --from=build /src/build/server/pricer_server /usr/local/bin/pricer_server
USER pricer
EXPOSE 8080
ENTRYPOINT ["pricer_server"]
CMD ["8080"]

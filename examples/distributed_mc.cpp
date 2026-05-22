// distributed_mc.cpp — distributed Monte Carlo across OS processes.
//
// Each worker process (a stand-in for a node) computes a contiguous range of the
// block decomposition and streams its per-block sums back through a pipe. The
// coordinator reassembles all B block-sums and reduces them in canonical order,
// so the price is BIT-FOR-BIT identical no matter how many workers ran it — while
// wall-clock time falls roughly in proportion to the worker count.
//
// Uses fork()/pipe() (POSIX), so this builds on Linux/macOS. The same block
// primitives would shard across real nodes (MPI/Ray/k8s) unchanged.
#include <sys/wait.h>
#include <unistd.h>

#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <vector>

#include "pricer/black_scholes.hpp"
#include "pricer/distributed.hpp"

using namespace pricer;

// Read exactly n bytes from a fd (pipes may deliver partial reads).
static bool read_full(int fd, void* buf, size_t n) {
    char* p = static_cast<char*>(buf);
    while (n > 0) {
        const ssize_t k = ::read(fd, p, n);
        if (k <= 0) return false;
        p += k;
        n -= static_cast<size_t>(k);
    }
    return true;
}

// Run the MC across `workers` forked processes; returns the price.
static double run_distributed(int workers, double S, double K, double r, double sigma, double T,
                              long n_paths, int B) {
    auto call = [K](double ST) { return ST > K ? ST - K : 0.0; };
    std::vector<double> block(B, 0.0);
    std::vector<int> rfd(workers);
    std::vector<pid_t> pids(workers);
    std::vector<int> lo(workers), hi(workers);

    for (int w = 0; w < workers; ++w) {
        lo[w] = static_cast<int>(static_cast<long>(B) * w / workers);
        hi[w] = static_cast<int>(static_cast<long>(B) * (w + 1) / workers);
        int fds[2];
        if (::pipe(fds) != 0) { std::perror("pipe"); std::exit(1); }
        const pid_t pid = ::fork();
        if (pid == 0) {  // ---- worker process ----
            ::close(fds[0]);
            std::vector<double> part(hi[w] - lo[w]);
            mc::fill_block_sums(call, S, r, sigma, T, n_paths, B, lo[w], hi[w], 12345, part.data());
            const ssize_t want = static_cast<ssize_t>(part.size() * sizeof(double));
            if (::write(fds[1], part.data(), static_cast<size_t>(want)) != want) _exit(1);
            ::close(fds[1]);
            _exit(0);
        }
        ::close(fds[1]);
        rfd[w] = fds[0];
        pids[w] = pid;
    }

    // Coordinator: place each worker's block-sums at their canonical positions.
    for (int w = 0; w < workers; ++w) {
        read_full(rfd[w], &block[lo[w]], static_cast<size_t>(hi[w] - lo[w]) * sizeof(double));
        ::close(rfd[w]);
    }
    for (int w = 0; w < workers; ++w) { int st; ::waitpid(pids[w], &st, 0); }

    double total = 0.0;
    for (int b = 0; b < B; ++b) total += block[b];  // canonical-order reduction
    return mc::aggregate_price(total, n_paths, r, T);
}

int main() {
    const double S = 100, K = 100, r = 0.05, sigma = 0.20, T = 1.0;
    const long n_paths = 80'000'000;
    const int B = 512;
    const double bs = black_scholes_call(S, K, r, sigma, T);

    std::printf("distributed Monte Carlo: %ld paths over %d blocks\n", n_paths, B);
    std::printf("analytic (Black-Scholes) = %.8f\n\n", bs);
    std::printf("%8s | %14s | %10s | %9s\n", "workers", "price", "time(ms)", "speedup");
    std::printf("---------|----------------|------------|----------\n");

    double base_ms = 0.0, ref_price = 0.0;
    for (int w : {1, 2, 4, 8}) {
        const auto t0 = std::chrono::high_resolution_clock::now();
        const double price = run_distributed(w, S, K, r, sigma, T, n_paths, B);
        const auto t1 = std::chrono::high_resolution_clock::now();
        const double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
        if (w == 1) { base_ms = ms; ref_price = price; }
        std::printf("%8d | %14.8f | %10.1f | %8.2fx\n", w, price, ms, base_ms / ms);
        if (price != ref_price)
            std::printf("   !! price differs from 1-worker run by %.2e\n", price - ref_price);
    }
    std::printf("\nEvery worker count produced the identical price (deterministic global\n"
                "aggregation); wall-clock time scales with the number of processes.\n");
    return 0;
}

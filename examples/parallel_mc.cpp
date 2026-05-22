// parallel_mc.cpp — speed up Monte Carlo across CPU cores.
// Runs the same total work single-threaded and across all cores, with a
// deterministic per-thread seed and a final reduction of the payoff sums.
#include "pricer/black_scholes.hpp"
#include <cstdio>
#include <random>
#include <chrono>
#include <thread>
#include <vector>

using namespace pricer;

// Returns the (undiscounted) payoff sum for `n` paths from a given seed.
static double payoff_sum(unsigned seed, long n, double S, double K,
                         double drift, double vol) {
    std::mt19937_64 rng(seed);
    std::normal_distribution<double> Z(0.0, 1.0);
    double sum = 0.0;
    for (long i = 0; i < n; ++i) {
        const double ST = S * std::exp(drift + vol * Z(rng));
        if (ST > K) sum += ST - K;
    }
    return sum;
}

int main() {
    const double S = 100, K = 100, r = 0.05, sigma = 0.20, T = 1.0;
    const long total = 200'000'000;
    const double drift = (r - 0.5 * sigma * sigma) * T;
    const double vol = sigma * std::sqrt(T);
    const double disc = std::exp(-r * T);
    const double exact = black_scholes_call(S, K, r, sigma, T);

    const unsigned ncores = std::thread::hardware_concurrency();
    std::printf("analytic = %.6f / total paths = %ld / logical cores = %u\n\n",
                exact, total, ncores);

    {  // single thread
        const auto t0 = std::chrono::high_resolution_clock::now();
        const double price = disc * payoff_sum(1, total, S, K, drift, vol) / total;
        const auto t1 = std::chrono::high_resolution_clock::now();
        const double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
        std::printf("1 thread   : price=%.6f  time=%8.1f ms\n", price, ms);
    }

    {  // all cores
        const auto t0 = std::chrono::high_resolution_clock::now();
        std::vector<std::thread> threads;
        std::vector<double> partial(ncores, 0.0);
        const long per = total / ncores;
        for (unsigned t = 0; t < ncores; ++t) {
            const long n = (t == ncores - 1) ? (total - per * (ncores - 1)) : per;
            threads.emplace_back([&, t, n] {
                partial[t] = payoff_sum(100 + t, n, S, K, drift, vol);
            });
        }
        for (auto& th : threads) th.join();
        double sum = 0.0;
        for (double p : partial) sum += p;
        const double price = disc * sum / total;
        const auto t1 = std::chrono::high_resolution_clock::now();
        const double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
        std::printf("%u threads  : price=%.6f  time=%8.1f ms\n", ncores, price, ms);
    }
    std::printf("\nSpeedup scales with core count (bounded by memory bandwidth etc.).\n");
    return 0;
}

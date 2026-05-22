// pricer/parallel.hpp — deterministic multithreaded Monte Carlo.
//
// `price_terminal_parallel` splits the work into a FIXED number of blocks (not a
// function of the thread count). Each block has a deterministic sub-seed and is
// summed independently; partial sums are combined in block order. The result is
// therefore bit-for-bit identical regardless of how many threads run it — a key
// property for reproducible risk numbers.
#pragma once
#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdint>
#include <random>
#include <thread>
#include <vector>

namespace pricer::mc {

namespace detail {
// SplitMix64: mixes a counter into a well-distributed 64-bit seed.
inline std::uint64_t splitmix64(std::uint64_t x) {
    x += 0x9E3779B97F4A7C15ull;
    x = (x ^ (x >> 30)) * 0xBF58476D1CE4E5B9ull;
    x = (x ^ (x >> 27)) * 0x94D049BB133111EBull;
    return x ^ (x >> 31);
}
}  // namespace detail

// Discounted mean of payoff(S_T) over `n_paths`, computed across `n_threads`
// (0 = hardware concurrency). Deterministic in `n_threads`.
template <class Payoff>
double price_terminal_parallel(Payoff payoff, double S, double r, double sigma, double T,
                               long n_paths, std::uint64_t seed = 12345, unsigned n_threads = 0) {
    if (n_threads == 0) n_threads = std::max(1u, std::thread::hardware_concurrency());
    constexpr int BLOCKS = 512;  // fixed work partition → result independent of thread count
    const double drift = (r - 0.5 * sigma * sigma) * T;
    const double vol = sigma * std::sqrt(T);
    const long base = n_paths / BLOCKS, rem = n_paths % BLOCKS;

    std::vector<double> partial(BLOCKS, 0.0);
    std::atomic<int> next{0};

    auto worker = [&] {
        int b;
        while ((b = next.fetch_add(1)) < BLOCKS) {
            const long start = static_cast<long>(b) * base + std::min<long>(b, rem);
            const long cnt = base + (b < rem ? 1 : 0);
            (void)start;  // block identity comes from its seed, not its index range
            std::mt19937_64 rng(detail::splitmix64(seed + static_cast<std::uint64_t>(b) + 1));
            std::normal_distribution<double> Z(0.0, 1.0);
            double s = 0.0;
            for (long i = 0; i < cnt; ++i) {
                const double ST = S * std::exp(drift + vol * Z(rng));
                s += payoff(ST);
            }
            partial[b] = s;
        }
    };

    std::vector<std::thread> threads;
    threads.reserve(n_threads);
    for (unsigned t = 0; t < n_threads; ++t) threads.emplace_back(worker);
    for (auto& th : threads) th.join();

    double sum = 0.0;
    for (int b = 0; b < BLOCKS; ++b) sum += partial[b];  // fixed-order reduction
    return std::exp(-r * T) * (sum / static_cast<double>(n_paths));
}

}  // namespace pricer::mc

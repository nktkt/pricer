// pricer/parallel_simd.hpp — multicore + SIMD Monte Carlo (both Phase-3 axes).
//
// Combines the deterministic block decomposition of parallel.hpp with the SIMD,
// counter-based path generation of simd_mc.hpp — the two single-node speedups,
// stacked. Path i still draws cb_normal(seed, i) for a GLOBAL counter i, so the
// price is independent of both the thread count and the SIMD lane width: the
// work partition fixes which counters land in which block, every block sums its
// own counters in a fixed order, and the blocks are reduced in index order, so
// the result is bit-identical no matter how many threads run it.
#pragma once
#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdint>
#include <thread>
#include <vector>

#include "pricer/rng.hpp"
#include "pricer/simd_mc.hpp"

namespace pricer::mc {

#if PRICER_HAVE_SIMD
namespace detail {

// Sum of payoff(S_T) over the global counter range [start, start+cnt): W paths
// per step via SIMD path generation, with a scalar tail for the final cnt mod W.
template <class Payoff>
double sum_range_cb_simd(Payoff payoff, double S, double drift, double vol, std::uint64_t seed,
                         long start, long cnt) {
    const simd::vd vS = simd::splat(S), vdrift = simd::splat(drift), vvol = simd::splat(vol);
    const long w = simd::kWidth;
    double s = 0.0;
    long i = 0;
    for (; i + w <= cnt; i += w) {
        simd::vd z = cb_normal_simd(seed, static_cast<std::uint64_t>(start + i));
        simd::vd ST = vS * simd::v_exp(vdrift + vvol * z);
        for (int j = 0; j < w; ++j) s += payoff(ST[j]);
    }
    for (; i < cnt; ++i)
        s += payoff(S * std::exp(drift + vol * cb_normal(seed,
                                                         static_cast<std::uint64_t>(start + i))));
    return s;
}

}  // namespace detail
#endif

// Discounted mean payoff over n_paths using multicore + SIMD path generation.
// Deterministic in n_threads (0 = hardware concurrency): the fixed block
// partition and fixed-order reduction make the result bit-identical for any
// thread count.
template <class Payoff>
double price_terminal_cb_parallel_simd(Payoff payoff, double S, double r, double sigma, double T,
                                       long n_paths, std::uint64_t seed = 12345,
                                       unsigned n_threads = 0) {
#if PRICER_HAVE_SIMD
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
            // Contiguous GLOBAL counter range for this block; every counter in
            // [0, n_paths) belongs to exactly one block, so draws are unchanged.
            const long start = static_cast<long>(b) * base + std::min<long>(b, rem);
            const long cnt = base + (b < rem ? 1 : 0);
            partial[b] = detail::sum_range_cb_simd(payoff, S, drift, vol, seed, start, cnt);
        }
    };

    std::vector<std::thread> threads;
    threads.reserve(n_threads);
    for (unsigned t = 0; t < n_threads; ++t) threads.emplace_back(worker);
    for (auto& th : threads) th.join();

    double sum = 0.0;
    for (int b = 0; b < BLOCKS; ++b) sum += partial[b];  // fixed-order reduction
    return std::exp(-r * T) * (sum / static_cast<double>(n_paths));
#else
    (void)n_threads;
    return price_terminal_cb_simd(payoff, S, r, sigma, T, n_paths, seed);  // serial scalar fallback
#endif
}

}  // namespace pricer::mc

// pricer/rng.hpp — counter-based RNG for SIMD / parallel / distributed paths.
//
// A counter-based generator has NO sequential state: the random draw for path i
// is a pure function of (seed, i). That makes path generation embarrassingly
// parallel and vectorization-friendly (a tight, branch-light loop over a batch
// of counters), and it makes results reproducible regardless of how the work is
// split across lanes, threads, or nodes — unlike a stateful engine such as
// std::mt19937_64.
#pragma once
#include <cmath>
#include <cstdint>

#include "pricer/parallel.hpp"  // mc::detail::splitmix64
#include "pricer/qmc.hpp"       // inv_norm_cdf

namespace pricer {

// Uniform in the open interval (0, 1) from (seed, counter). Pure integer mixing,
// so it vectorizes and is identical on every platform.
inline double cb_uniform(std::uint64_t seed, std::uint64_t counter) {
    const std::uint64_t z = mc::detail::splitmix64(seed + counter * 0x9E3779B97F4A7C15ull);
    return ((z >> 11) + 0.5) * (1.0 / 9007199254740992.0);  // 53-bit mantissa, never 0 or 1
}

// Standard normal from (seed, counter) via the inverse CDF.
inline double cb_normal(std::uint64_t seed, std::uint64_t counter) {
    return inv_norm_cdf(cb_uniform(seed, counter));
}

namespace mc {

// Fill `out[0..W)` with terminal prices for counters [base, base+W) under GBM.
// A tight, branch-light loop the compiler can vectorize.
inline void cb_terminal_batch(double S, double drift, double vol, std::uint64_t seed,
                              std::uint64_t base, int W, double* out) {
    for (int j = 0; j < W; ++j) out[j] = S * std::exp(drift + vol * cb_normal(seed, base + j));
}

// Terminal-value MC using the counter-based RNG. Reproducible for a given n and
// seed regardless of how the counters are partitioned.
template <class Payoff>
double price_terminal_cb(Payoff payoff, double S, double r, double sigma, double T, long n,
                         std::uint64_t seed = 12345) {
    const double drift = (r - 0.5 * sigma * sigma) * T;
    const double vol = sigma * std::sqrt(T);
    double sum = 0.0;
    for (long i = 0; i < n; ++i)
        sum += payoff(S * std::exp(drift + vol * cb_normal(seed, static_cast<std::uint64_t>(i))));
    return std::exp(-r * T) * (sum / static_cast<double>(n));
}

}  // namespace mc
}  // namespace pricer

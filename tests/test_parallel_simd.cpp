// Tests for the combined multicore + SIMD Monte Carlo engine (parallel_simd.hpp).
#include "check.hpp"
#include "pricer/black_scholes.hpp"
#include "pricer/parallel_simd.hpp"
#include "pricer/simd_mc.hpp"

using namespace pricer;

int main() {
    const double S = 100, K = 100, r = 0.05, sigma = 0.20, T = 1.0;
    auto call = [K](double ST) { return ST > K ? ST - K : 0.0; };
    const long n = 4'000'000;
    const double bs = black_scholes_call(S, K, r, sigma, T);

    // Converges to the analytic price.
    const double p8 = mc::price_terminal_cb_parallel_simd(call, S, r, sigma, T, n, 12345, 8);
    check::approx("parallel+SIMD MC vs Black-Scholes", p8, bs, 0.02);

#if PRICER_HAVE_SIMD
    // Deterministic in the thread count: bit-identical for 1, 2, 4, 8 threads.
    const double p1 = mc::price_terminal_cb_parallel_simd(call, S, r, sigma, T, n, 12345, 1);
    const double p2 = mc::price_terminal_cb_parallel_simd(call, S, r, sigma, T, n, 12345, 2);
    const double p4 = mc::price_terminal_cb_parallel_simd(call, S, r, sigma, T, n, 12345, 4);
    check::is_true("thread count 1==2 (bit-identical)", p1 == p2);
    check::is_true("thread count 1==4 (bit-identical)", p1 == p4);
    check::is_true("thread count 1==8 (bit-identical)", p1 == p8);

    // Auto thread count (0 = hardware concurrency) gives the same price.
    const double p0 = mc::price_terminal_cb_parallel_simd(call, S, r, sigma, T, n, 12345, 0);
    check::is_true("auto threads == explicit", p0 == p1);

    // Agrees with the serial SIMD engine to tolerance (block boundaries reshuffle
    // the SIMD grouping, so it is not bit-identical, only floating-point close).
    const double serial = mc::price_terminal_cb_simd(call, S, r, sigma, T, n);
    check::approx("parallel+SIMD vs serial SIMD", p1, serial, 1e-9);

    // A path count not divisible by the lane width still prices correctly.
    const double odd = mc::price_terminal_cb_parallel_simd(call, S, r, sigma, T, n + 7, 12345, 4);
    check::approx("parallel+SIMD handles ragged n", odd, bs, 0.02);
#else
    check::is_true("SIMD unavailable: parallel+SIMD falls back", true);
#endif
    return check::report("parallel_simd");
}

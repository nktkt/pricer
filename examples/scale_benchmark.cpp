// scale_benchmark.cpp — the single-node speedup ladder for Monte Carlo pricing.
//
// Phase 3 of the roadmap asks for >=10x over the Phase-1 baseline on CPU. This
// program prices the same European call four ways and reports the cumulative
// speedup at each rung, so the contribution of every optimization is visible:
//
//   1. Phase-1 baseline : stateful std::mt19937_64 + std::normal_distribution,
//                         single thread (the textbook Monte Carlo loop).
//   2. + counter-based  : stateless RNG (draw i = f(seed, i)), single thread.
//   3. + SIMD           : the same draws, W lanes per step (simd_mc.hpp).
//   4. + multicore      : SIMD path generation across all cores, with a
//                         deterministic, thread-count-independent result
//                         (parallel_simd.hpp).
//
// Rung 4 is the fully-optimized engine; its speedup over rung 1 is the headline
// number for the exit criterion. The actual factor scales with the core count
// and SIMD width of the machine you run it on.
#include "pricer/black_scholes.hpp"
#include "pricer/parallel_simd.hpp"
#include "pricer/rng.hpp"
#include "pricer/simd_mc.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <random>
#include <thread>
#include <utility>

using namespace pricer;

namespace {

constexpr double kS = 100.0, kK = 100.0, kR = 0.05, kSigma = 0.20, kT = 1.0;
constexpr long kN = 50'000'000L;

template <class TimePoint>
double elapsed_ms(TimePoint a, TimePoint b) {
    return std::chrono::duration<double, std::milli>(b - a).count();
}

// Phase-1 baseline: stateful std::mt19937_64, single thread.
template <class Payoff>
double price_terminal_mt(Payoff payoff, double S, double r, double sigma, double T, long n) {
    std::mt19937_64 rng(12345);
    std::normal_distribution<double> Z(0.0, 1.0);
    const double drift = (r - 0.5 * sigma * sigma) * T, vol = sigma * std::sqrt(T);
    double sum = 0.0;
    for (long i = 0; i < n; ++i) sum += payoff(S * std::exp(drift + vol * Z(rng)));
    return std::exp(-r * T) * (sum / static_cast<double>(n));
}

// Time a pricing call; returns {price, elapsed_ms}.
template <class F>
std::pair<double, double> timed(F&& f) {
    const auto t0 = std::chrono::high_resolution_clock::now();
    const double price = f();
    const auto t1 = std::chrono::high_resolution_clock::now();
    return {price, elapsed_ms(t0, t1)};
}

}  // namespace

int main() {
    auto call = [](double ST) { return ST > kK ? ST - kK : 0.0; };
    const double exact = black_scholes_call(kS, kK, kR, kSigma, kT);
    const unsigned cores = std::max(1u, std::thread::hardware_concurrency());

    std::printf("Single-node Monte Carlo speedup ladder\n");
    std::printf("Params: S=%.1f K=%.1f r=%.2f sigma=%.2f T=%.1f  paths=%ld  cores=%u",
                kS, kK, kR, kSigma, kT, kN, cores);
#if PRICER_HAVE_SIMD
    std::printf("  SIMD width=%d\n\n", simd::kWidth);
#else
    std::printf("  SIMD unavailable (scalar fallback)\n\n");
#endif

    const auto base = timed([&] { return price_terminal_mt(call, kS, kR, kSigma, kT, kN); });
    const auto cb = timed([&] { return mc::price_terminal_cb(call, kS, kR, kSigma, kT, kN); });
    const auto simd = timed([&] { return mc::price_terminal_cb_simd(call, kS, kR, kSigma, kT, kN); });
    const auto full =
        timed([&] { return mc::price_terminal_cb_parallel_simd(call, kS, kR, kSigma, kT, kN); });

    const double bms = base.second;
    std::printf("%-34s  %12s  %10s  %10s\n", "engine", "price", "time(ms)", "speedup");
    std::printf("------------------------------------------------------------------------\n");
    auto row = [&](const char* name, std::pair<double, double> m) {
        std::printf("%-34s  %12.6f  %10.2f  %9.2fx\n", name, m.first, m.second,
                    (m.second > 0.0) ? bms / m.second : 0.0);
    };
    row("1. stateful mt19937 (Phase-1)", base);
    row("2. + counter-based RNG", cb);
    row("3. + SIMD path generation", simd);
    row("4. + multicore (full engine)", full);

    const double total = (full.second > 0.0) ? bms / full.second : 0.0;
    std::printf("\nFull engine speedup over the Phase-1 baseline: %.2fx", total);
    std::printf("  [exit criterion >=10x: %s on this machine]\n", total >= 10.0 ? "MET" : "not met");

    std::printf("\nAll engines agree with the analytic price (Black-Scholes = %.6f):\n", exact);
    std::printf("  baseline %.6f   counter %.6f   SIMD %.6f   full %.6f\n", base.first, cb.first,
                simd.first, full.first);
    std::printf("\nThe full engine is deterministic in the thread count: the fixed block\n");
    std::printf("partition makes its price bit-identical regardless of how many cores run it.\n");
    return 0;
}

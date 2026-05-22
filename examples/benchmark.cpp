// benchmark.cpp — Monte Carlo performance harness (throughput baseline).
//
// This program is NOT about pricing accuracy (the unit tests cover that); it
// exists to measure how many GBM paths per second the engine can simulate, so
// that future optimization work (Phase 3 of the roadmap) can be tracked against
// a stable, reproducible baseline.
//
// Two sections are reported:
//   1. Single-threaded throughput across a range of path counts. This drives
//      `pricer::mc::price_terminal` directly, so it benchmarks the shipped API.
//   2. Multi-thread scaling for a fixed large workload, split across 1/2/4/N
//      threads, reporting throughput and speedup vs. the single-thread case.
//      The threaded path uses a small local helper that mirrors the engine's
//      GBM terminal simulation, because `price_terminal` is single-threaded by
//      design and we want each worker to own an independent RNG stream.
//
// Throughput is reported in millions of paths per second (Mpaths/s), which is
// the most meaningful unit for comparing runs over time.

#include "pricer/black_scholes.hpp"
#include "pricer/monte_carlo.hpp"

#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <random>
#include <thread>
#include <vector>

using namespace pricer;

// ---------------------------------------------------------------------------
// Shared market parameters and the instrument under test.
// ---------------------------------------------------------------------------
namespace {

constexpr double kS     = 100.0;  // spot
constexpr double kK     = 100.0;  // strike
constexpr double kR     = 0.05;   // risk-free rate
constexpr double kSigma = 0.20;   // volatility
constexpr double kT     = 1.0;    // time to expiry (years)

// European call payoff: max(S_T - K, 0). Captured by value where needed.
inline double call_payoff(double ST) { return ST > kK ? ST - kK : 0.0; }

// Convenience: elapsed time in milliseconds between two clock samples.
template <class TimePoint>
inline double elapsed_ms(TimePoint t0, TimePoint t1) {
    return std::chrono::duration<double, std::milli>(t1 - t0).count();
}

// Throughput in Mpaths/s given a path count and an elapsed time in ms.
inline double throughput_mpps(long n_paths, double ms) {
    if (ms <= 0.0) return 0.0;
    // paths / second = n_paths / (ms / 1000); divide by 1e6 for "mega".
    return static_cast<double>(n_paths) / (ms / 1000.0) / 1.0e6;
}

// ---------------------------------------------------------------------------
// Local GBM payoff-sum kernel for the threaded scaling test.
//
// Returns the *undiscounted* sum of payoffs over `n` paths, each drawn from an
// independent N(0,1) stream seeded by `seed`. Discounting and averaging are
// done once, after all partial sums are reduced, so the math matches
// `price_terminal` exactly (sum -> mean -> discount).
// ---------------------------------------------------------------------------
double payoff_sum(std::uint64_t seed, long n, double drift, double vol) {
    std::mt19937_64 rng(seed);
    std::normal_distribution<double> Z(0.0, 1.0);
    double sum = 0.0;
    for (long i = 0; i < n; ++i) {
        // S_T = S * exp((r - 0.5*sigma^2)*T + sigma*sqrt(T)*Z)
        const double ST = kS * std::exp(drift + vol * Z(rng));
        sum += call_payoff(ST);
    }
    return sum;
}

// ---------------------------------------------------------------------------
// Section 1 — single-threaded throughput across several path counts.
// ---------------------------------------------------------------------------
void run_single_threaded_section(double exact) {
    std::printf("== Single-threaded throughput (pricer::mc::price_terminal) ==\n");
    std::printf("%12s  %10s  %14s  %12s  %12s\n",
                "paths", "time(ms)", "Mpaths/s", "price", "abs err");
    std::printf("------------------------------------------------------------------------\n");

    // A few representative workload sizes. Larger sizes amortize the per-call
    // setup (RNG construction, drift/vol precompute) and give a steadier number.
    const long path_counts[] = {1'000'000L, 10'000'000L, 50'000'000L};

    for (long n : path_counts) {
        // One warm-up at a small fixed size keeps caches/branch predictors warm
        // without skewing the measured run for this particular size.
        volatile double warm = mc::price_terminal(call_payoff, kS, kR, kSigma, kT,
                                                  100'000L, /*seed=*/7);
        (void)warm;

        const auto t0 = std::chrono::high_resolution_clock::now();
        const double price = mc::price_terminal(call_payoff, kS, kR, kSigma, kT, n);
        const auto t1 = std::chrono::high_resolution_clock::now();

        const double ms  = elapsed_ms(t0, t1);
        const double mpps = throughput_mpps(n, ms);
        const double err = std::fabs(price - exact);

        std::printf("%12ld  %10.2f  %14.2f  %12.6f  %12.6f\n",
                    n, ms, mpps, price, err);
    }
    std::printf("\n");
}

// ---------------------------------------------------------------------------
// Section 2 — multi-thread scaling for a fixed large total path count.
// ---------------------------------------------------------------------------
void run_scaling_section(double exact) {
    const long total = 100'000'000L;  // fixed total work shared by all threads
    const double drift = (kR - 0.5 * kSigma * kSigma) * kT;
    const double vol   = kSigma * std::sqrt(kT);
    const double disc  = std::exp(-kR * kT);

    std::printf("== Multi-thread scaling (fixed total = %ld paths) ==\n", total);
    std::printf("%9s  %10s  %14s  %10s  %12s  %12s\n",
                "threads", "time(ms)", "Mpaths/s", "speedup", "price", "abs err");
    std::printf("------------------------------------------------------------------------------\n");

    // Candidate thread counts: 1, 2, 4, and the machine's logical core count.
    const unsigned hw = std::thread::hardware_concurrency();
    unsigned thread_counts[] = {1u, 2u, 4u, hw == 0 ? 1u : hw};

    double baseline_ms = 0.0;  // time for the 1-thread run, used for speedup

    for (unsigned nt : thread_counts) {
        if (nt == 0) nt = 1;  // defensive: hardware_concurrency may report 0

        const auto t0 = std::chrono::high_resolution_clock::now();

        // Split `total` into nt chunks; the last chunk absorbs any remainder.
        std::vector<std::thread> workers;
        std::vector<double> partial(nt, 0.0);
        const long per = total / static_cast<long>(nt);

        for (unsigned t = 0; t < nt; ++t) {
            const long n = (t == nt - 1)
                               ? (total - per * static_cast<long>(nt - 1))
                               : per;
            // Distinct seed per worker so the N(0,1) streams do not overlap.
            const std::uint64_t seed = 1000ull + t;
            workers.emplace_back([&partial, t, seed, n, drift, vol] {
                partial[t] = payoff_sum(seed, n, drift, vol);
            });
        }
        for (auto& w : workers) w.join();

        // Reduce partial sums, then average and discount exactly once.
        double sum = 0.0;
        for (double p : partial) sum += p;
        const double price = disc * (sum / static_cast<double>(total));

        const auto t1 = std::chrono::high_resolution_clock::now();
        const double ms = elapsed_ms(t0, t1);

        if (nt == 1) baseline_ms = ms;  // first row establishes the baseline
        const double speedup = (baseline_ms > 0.0) ? baseline_ms / ms : 1.0;
        const double mpps = throughput_mpps(total, ms);
        const double err  = std::fabs(price - exact);

        std::printf("%9u  %10.2f  %14.2f  %10.2f  %12.6f  %12.6f\n",
                    nt, ms, mpps, speedup, price, err);
    }
    std::printf("\n");
}

}  // namespace

int main() {
    // Analytic reference price; every Monte Carlo estimate is compared to this.
    const double exact = black_scholes_call(kS, kK, kR, kSigma, kT);

    std::printf("Pricer Monte Carlo benchmark\n");
    std::printf("Params: S=%.1f K=%.1f r=%.2f sigma=%.2f T=%.1f\n",
                kS, kK, kR, kSigma, kT);
    std::printf("Analytic (Black-Scholes call) = %.6f\n", exact);
    std::printf("Logical cores reported        = %u\n\n",
                std::thread::hardware_concurrency());

    run_single_threaded_section(exact);
    run_scaling_section(exact);

    std::printf("Note: numbers are a baseline for tracking Phase 3 optimizations.\n");
    return 0;
}

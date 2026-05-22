// black_scholes_demo.cpp — price a European call two ways and compare.
//   (1) Black–Scholes closed form (analytic)
//   (2) Monte Carlo via the generic engine, with the payoff as a lambda
// Built on the `pricer` library headers.
#include "pricer/black_scholes.hpp"
#include "pricer/monte_carlo.hpp"
#include <cstdio>
#include <chrono>

using namespace pricer;

int main() {
    const double S = 100, K = 100, r = 0.05, sigma = 0.20, T = 1.0;
    const long n_paths = 10'000'000;

    const double bs = black_scholes_call(S, K, r, sigma, T);

    const auto t0 = std::chrono::high_resolution_clock::now();
    const double mc = mc::price_terminal(
        [K](double ST) { return ST > K ? ST - K : 0.0; }, S, r, sigma, T, n_paths);
    const auto t1 = std::chrono::high_resolution_clock::now();
    const double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();

    std::printf("S=%.1f K=%.1f r=%.2f sigma=%.2f T=%.1f  paths=%ld\n\n",
                S, K, r, sigma, T, n_paths);
    std::printf("(1) Black-Scholes : %.6f\n", bs);
    std::printf("(2) Monte Carlo   : %.6f\n", mc);
    std::printf("    error         : %.6f (%.4f%%)\n", mc - bs, (mc - bs) / bs * 100.0);
    std::printf("    MC time        : %.1f ms\n", ms);
    return 0;
}

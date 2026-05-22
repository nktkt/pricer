// convergence.cpp — accuracy vs. speed.
// Grow the number of Monte Carlo paths and watch the error against the analytic
// price shrink like 1/sqrt(N) (×100 paths → ~1/10 error).
#include "pricer/black_scholes.hpp"
#include <cstdio>
#include <random>
#include <chrono>

using namespace pricer;

static double monte_carlo_call(double S, double K, double r, double sigma, double T,
                               long n_paths) {
    std::mt19937_64 rng(12345);
    std::normal_distribution<double> Z(0.0, 1.0);
    const double drift = (r - 0.5 * sigma * sigma) * T;
    const double vol = sigma * std::sqrt(T);
    double sum = 0.0;
    for (long i = 0; i < n_paths; ++i) {
        const double ST = S * std::exp(drift + vol * Z(rng));
        if (ST > K) sum += ST - K;
    }
    return std::exp(-r * T) * (sum / static_cast<double>(n_paths));
}

int main() {
    const double S = 100, K = 100, r = 0.05, sigma = 0.20, T = 1.0;
    const double exact = black_scholes_call(S, K, r, sigma, T);
    std::printf("analytic (Black-Scholes) = %.6f\n\n", exact);
    std::printf("%12s | %12s | %10s | %10s\n", "paths", "MC price", "error", "ms");
    std::printf("-------------|--------------|------------|----------\n");

    for (long n = 10'000; n <= 100'000'000; n *= 10) {
        const auto t0 = std::chrono::high_resolution_clock::now();
        const double mc = monte_carlo_call(S, K, r, sigma, T, n);
        const auto t1 = std::chrono::high_resolution_clock::now();
        const double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
        std::printf("%12ld | %12.6f | %10.6f | %10.1f\n", n, mc, mc - exact, ms);
    }
    std::printf("\nEach ×100 in paths cuts the error to roughly 1/10 (the 1/sqrt(N) law).\n");
    return 0;
}

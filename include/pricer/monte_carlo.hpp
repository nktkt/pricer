// pricer/monte_carlo.hpp — generic Monte Carlo pricing under geometric Brownian motion.
#pragma once
#include <cmath>
#include <cstdint>
#include <random>

namespace pricer::mc {

// Price a European-style instrument by simulating the terminal spot
//   S_T = S * exp((r - 0.5*sigma^2)*T + sigma*sqrt(T)*Z),   Z ~ N(0,1)
// applying `payoff(S_T)` to each path, and discounting the average.
//
// `payoff` is any callable `double(double S_T)`, so the caller decides the
// instrument (call, put, digital, …) without this engine knowing about it.
template <class Payoff>
double price_terminal(Payoff&& payoff, double S, double r, double sigma, double T,
                      long n_paths, std::uint64_t seed = 12345) {
    std::mt19937_64 rng(seed);
    std::normal_distribution<double> Z(0.0, 1.0);
    const double drift = (r - 0.5 * sigma * sigma) * T;
    const double vol = sigma * std::sqrt(T);

    double sum = 0.0;
    for (long i = 0; i < n_paths; ++i) {
        const double ST = S * std::exp(drift + vol * Z(rng));
        sum += payoff(ST);
    }
    return std::exp(-r * T) * (sum / static_cast<double>(n_paths));
}

}  // namespace pricer::mc

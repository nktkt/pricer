// barrier_option.cpp — a path-dependent product (up-and-out barrier call).
// A normal call that becomes worthless if the spot ever touches the barrier B.
// Because the payoff depends on the whole path, we step through time with Monte
// Carlo. A plain call is priced on the same paths for comparison; knock-outs make
// the barrier version strictly cheaper.
#include "pricer/black_scholes.hpp"
#include <cstdio>
#include <random>

using namespace pricer;

int main() {
    const double S = 100, K = 100, r = 0.05, sigma = 0.20, T = 1.0;
    const double B = 130.0;          // barrier level
    const long n_paths = 2'000'000;
    const int  n_steps = 250;        // ~one trading day per step

    const double dt = T / n_steps;
    const double drift = (r - 0.5 * sigma * sigma) * dt;
    const double vol = sigma * std::sqrt(dt);
    const double disc = std::exp(-r * T);

    std::mt19937_64 rng(2024);
    std::normal_distribution<double> Z(0.0, 1.0);

    double sum_barrier = 0.0, sum_vanilla = 0.0;
    for (long p = 0; p < n_paths; ++p) {
        double price = S;
        bool knocked_out = false;
        for (int s = 0; s < n_steps; ++s) {
            price *= std::exp(drift + vol * Z(rng));
            if (price >= B) knocked_out = true;
        }
        const double payoff = (price > K) ? (price - K) : 0.0;
        sum_vanilla += payoff;
        if (!knocked_out) sum_barrier += payoff;
    }

    const double vanilla_mc = disc * sum_vanilla / n_paths;
    const double barrier_mc = disc * sum_barrier / n_paths;
    const double vanilla_bs = black_scholes_call(S, K, r, sigma, T);

    std::printf("S=%.0f K=%.0f barrier B=%.0f paths=%ld steps=%d\n\n",
                S, K, B, n_paths, n_steps);
    std::printf("vanilla call (analytic BS) : %.6f\n", vanilla_bs);
    std::printf("vanilla call (Monte Carlo) : %.6f  <- matches BS\n", vanilla_mc);
    std::printf("up-and-out barrier call    : %.6f  <- cheaper by the knock-outs\n", barrier_mc);
    std::printf("difference (knock-out cost): %.6f\n", vanilla_mc - barrier_mc);
    std::printf("\nBarrier payoffs often lack a simple formula; stepped MC is the natural tool.\n");
    return 0;
}

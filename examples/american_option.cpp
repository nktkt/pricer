// american_option.cpp — pricing an American option two ways, and the
// early-exercise premium.
//
// An American option can be exercised any time before expiry, so it is worth at
// least its European counterpart; the difference is the early-exercise premium.
// This demo prices an American put with a Cox–Ross–Rubinstein binomial tree and
// with Longstaff–Schwartz least-squares Monte Carlo, and cross-checks both
// against the European value (Black–Scholes / the binomial European mode). For a
// non-dividend stock the American *call* equals the European call, which the demo
// also shows.
#include "pricer/american.hpp"
#include "pricer/black_scholes.hpp"

#include <cstdio>

using namespace pricer;

int main() {
    const double S = 100, K = 100, r = 0.05, sigma = 0.20, T = 1.0;
    const int steps = 1000;
    const long n_paths = 200'000;
    const int n_steps = 50;

    std::printf("American option pricing   S=%.0f K=%.0f r=%.2f sigma=%.2f T=%.1f\n\n",
                S, K, r, sigma, T);

    // --- Put: European vs American, by tree and by LSM Monte Carlo ---
    const double eur_put = black_scholes_put(S, K, r, sigma, T);
    const double bin_put = binomial_price(OptionType::Put, S, K, r, sigma, T, steps, true);
    const double lsm_put = lsm_american(OptionType::Put, S, K, r, sigma, T, n_paths, n_steps);

    std::printf("PUT\n");
    std::printf("  European (Black-Scholes)        : %.4f\n", eur_put);
    std::printf("  American (binomial, %d steps)  : %.4f\n", steps, bin_put);
    std::printf("  American (LSM MC, %ld paths)  : %.4f\n", n_paths, lsm_put);
    std::printf("  early-exercise premium          : %.4f  (tree)\n", bin_put - eur_put);
    std::printf("  LSM vs binomial                 : %+.4f\n\n", lsm_put - bin_put);

    // --- Call on a non-dividend stock: American == European ---
    const double eur_call = black_scholes_call(S, K, r, sigma, T);
    const double bin_call = binomial_price(OptionType::Call, S, K, r, sigma, T, steps, true);
    std::printf("CALL (no dividends — early exercise never optimal)\n");
    std::printf("  European (Black-Scholes)        : %.4f\n", eur_call);
    std::printf("  American (binomial, %d steps)  : %.4f\n", steps, bin_call);
    std::printf("  difference                      : %.2e\n\n", bin_call - eur_call);

    // --- Deep in-the-money American put ~ intrinsic (exercise immediately) ---
    const double deep = binomial_price(OptionType::Put, 60.0, K, r, sigma, T, steps, true);
    std::printf("Deep-ITM American put (S=60, K=100): %.4f   (intrinsic K-S = %.1f)\n", deep, 40.0);
    return 0;
}

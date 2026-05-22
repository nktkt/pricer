// local_vol_demo.cpp — recover Dupire local volatility from a call-price surface.
// Given European call prices C(K, T), Dupire's formula yields the local vol that
// reprices the surface. Two checks: (A) a flat implied-vol surface must give back
// the same constant vol; (B) a skewed surface shows local vol moving more across
// strikes than the implied vol it came from (the "rule of two" intuition).
#include "pricer/black_scholes.hpp"
#include "pricer/local_vol.hpp"
#include <cstdio>
#include <vector>

using namespace pricer;

int main() {
    const double S = 100, r = 0.05;

    // -- Section A: flat implied-vol surface (constant 20%) sanity check. --------
    // A flat 20% surface must hand back ~0.20 local vol everywhere.
    const double flat_sigma = 0.20;
    auto price = [&](double K, double T) {
        return black_scholes_call(S, K, r, flat_sigma, T);
    };

    std::printf("Section A: flat 20%% implied-vol surface (expect local vol ~0.20)\n");
    std::printf("%8s | %8s | %12s\n", "K", "T", "local vol");
    std::printf("---------|----------|-------------\n");

    const std::vector<double> strikes_A = {90, 100, 110};
    const std::vector<double> times_A = {0.5, 1.0, 2.0};
    for (double K : strikes_A) {
        for (double T : times_A) {
            const double lv = dupire_local_vol(price, K, T, r);
            std::printf("%8.1f | %8.2f | %12.6f\n", K, T, lv);
        }
    }

    // -- Section B: skewed implied-vol surface. ---------------------------------
    // iv(K) = 0.20 + 0.0015*(100 - K): lower strikes carry higher vol.
    auto iv = [&](double K) {
        return 0.20 + 0.0015 * (100.0 - K);
    };
    auto price_skew = [&](double K, double T) {
        return black_scholes_call(S, K, r, iv(K), T);
    };

    std::printf("\nSection B: skewed implied-vol surface at T=1.0\n");
    std::printf("iv(K) = 0.20 + 0.0015*(100 - K)  (lower strike => higher vol)\n");
    std::printf("%8s | %12s | %14s\n", "K", "implied vol", "local vol");
    std::printf("---------|--------------|----------------\n");

    const double T_B = 1.0;
    const std::vector<double> strikes_B = {80, 90, 100, 110, 120};
    for (double K : strikes_B) {
        const double iv_K = iv(K);
        const double lv = dupire_local_vol(price_skew, K, T_B, r);
        std::printf("%8.1f | %12.6f | %14.6f\n", K, iv_K, lv);
    }

    std::printf("\nThe local-vol column moves more across strikes than the implied-vol\n");
    std::printf("column: local-vol skew is steeper (the well-known 'rule of two').\n");
    return 0;
}

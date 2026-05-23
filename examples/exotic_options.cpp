// exotic_options.cpp — pricing path-dependent exotics (Asian, barrier, lookback)
// two ways each: a closed form and a Monte Carlo engine that should agree.
//
// Asian options pay on the *average* spot, barriers knock in/out when the spot
// crosses a level, and lookbacks pay on the path's running extreme. The closed
// forms are exact (geometric Asian) or exact under continuous monitoring
// (barrier, lookback); the Monte Carlo engine steps the path with the
// counter-based RNG and, for barriers/lookbacks, applies the
// Broadie–Glasserman–Kou continuity correction so discrete stepping converges to
// the continuous price.
#include "pricer/black_scholes.hpp"
#include "pricer/exotics.hpp"

#include <cstdio>

using namespace pricer;

int main() {
    const double S = 100, K = 100, r = 0.05, sigma = 0.20, T = 1.0;
    const int n_steps = 250;
    const long n_paths = 500'000;

    std::printf("Exotic options   S=%.0f K=%.0f r=%.2f sigma=%.2f T=%.1f   (%ld paths, %d steps)\n\n",
                S, K, r, sigma, T, n_paths, n_steps);

    // --- Asian (average price, fixed strike) ---
    const double geo_cf = geometric_asian_price(OptionType::Call, S, K, r, sigma, T, n_steps);
    const double geo_mc =
        asian_price_mc(OptionType::Call, AverageType::Geometric, S, K, r, sigma, T, n_steps, n_paths);
    const double ari_mc =
        asian_price_mc(OptionType::Call, AverageType::Arithmetic, S, K, r, sigma, T, n_steps, n_paths);
    std::printf("ASIAN call (average price)\n");
    std::printf("  geometric  (closed form)   : %.4f\n", geo_cf);
    std::printf("  geometric  (Monte Carlo)   : %.4f   <- matches closed form\n", geo_mc);
    std::printf("  arithmetic (Monte Carlo)   : %.4f   <- dearer (AM >= GM)\n", ari_mc);
    std::printf("  vanilla European (BS)      : %.4f   <- averaging dampens vol\n\n",
                black_scholes_call(S, K, r, sigma, T));

    // --- Barrier (single barrier, continuous monitoring) ---
    const double B = 130.0;
    const double uo_cf = barrier_price(OptionType::Call, BarrierType::UpOut, S, K, B, r, sigma, T);
    const double ui_cf = barrier_price(OptionType::Call, BarrierType::UpIn, S, K, B, r, sigma, T);
    const double uo_mc =
        barrier_price_mc(OptionType::Call, BarrierType::UpOut, S, K, B, r, sigma, T, n_steps, n_paths);
    std::printf("BARRIER call (up barrier B=%.0f)\n", B);
    std::printf("  up-and-out (closed form)   : %.4f\n", uo_cf);
    std::printf("  up-and-out (Monte Carlo)   : %.4f   <- BGK-corrected, matches\n", uo_mc);
    std::printf("  up-and-in  (closed form)   : %.4f\n", ui_cf);
    std::printf("  in + out                   : %.4f   <- equals vanilla %.4f (parity)\n\n",
                uo_cf + ui_cf, black_scholes_call(S, K, r, sigma, T));

    // --- Lookback (floating strike, freshly issued) ---
    const double lc_cf = lookback_floating_price(OptionType::Call, S, r, sigma, T);
    const double lp_cf = lookback_floating_price(OptionType::Put, S, r, sigma, T);
    const double lc_mc = lookback_floating_price_mc(OptionType::Call, S, r, sigma, T, n_steps, n_paths);
    const double lp_mc = lookback_floating_price_mc(OptionType::Put, S, r, sigma, T, n_steps, n_paths);
    std::printf("LOOKBACK (floating strike)\n");
    std::printf("  call S_T - min  (closed)   : %.4f\n", lc_cf);
    std::printf("  call            (MC)       : %.4f\n", lc_mc);
    std::printf("  put  max - S_T  (closed)   : %.4f\n", lp_cf);
    std::printf("  put             (MC)       : %.4f\n", lp_mc);
    return 0;
}

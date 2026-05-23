// pricer/digital.hpp — digital (binary) options.
//
// A digital option pays a discontinuous amount that depends only on whether the
// option finishes in the money:
//   * cash-or-nothing pays a fixed amount `cash` if in the money, else nothing;
//   * asset-or-nothing pays the underlying S_T itself if in the money.
// Both have simple closed forms (the discounted in-the-money probability times
// the payoff), and together they decompose the vanilla option exactly:
//   vanilla call = asset-or-nothing call − K · cash-or-nothing call (cash = 1),
// which is the cleanest cross-check for the formulas below. A counter-based Monte
// Carlo engine prices the same discontinuous payoffs for comparison.
//
// All prices carry a continuous dividend yield q (q = 0 by default).
#pragma once
#include <cmath>
#include <cstdint>

#include "pricer/black_scholes.hpp"  // OptionType, norm_cdf (via normal.hpp)
#include "pricer/rng.hpp"            // cb_normal

namespace pricer {

enum class DigitalType { CashOrNothing, AssetOrNothing };

// Closed-form price of a digital option. A cash-or-nothing option pays `cash`
// when in the money (value = cash · e^{-rT} · N(±d2), the discounted exercise
// probability); an asset-or-nothing option pays S_T (value = S · e^{-qT} · N(±d1)).
inline double digital_price(OptionType type, DigitalType kind, double S, double K, double r,
                            double sigma, double T, double cash = 1.0, double q = 0.0) {
    const double sqrtT = std::sqrt(T);
    const double d1 = (std::log(S / K) + (r - q + 0.5 * sigma * sigma) * T) / (sigma * sqrtT);
    const double d2 = d1 - sigma * sqrtT;
    const double disc = std::exp(-r * T);
    const double discq = std::exp(-q * T);
    if (kind == DigitalType::CashOrNothing)
        return type == OptionType::Call ? cash * disc * norm_cdf(d2) : cash * disc * norm_cdf(-d2);
    return type == OptionType::Call ? S * discq * norm_cdf(d1) : S * discq * norm_cdf(-d1);
}

// Monte Carlo price of a digital option, simulating the terminal spot with the
// counter-based RNG (reproducible). The payoff is discontinuous at the strike,
// so the estimator has higher variance than a vanilla one.
inline double digital_price_mc(OptionType type, DigitalType kind, double S, double K, double r,
                               double sigma, double T, long n_paths, std::uint64_t seed = 12345,
                               double cash = 1.0, double q = 0.0) {
    if (n_paths < 1) n_paths = 1;
    const double drift = (r - q - 0.5 * sigma * sigma) * T;
    const double vol = sigma * std::sqrt(T);
    const double disc = std::exp(-r * T);
    double sum = 0.0;
    for (long p = 0; p < n_paths; ++p) {
        const double ST = S * std::exp(drift + vol * cb_normal(seed, static_cast<std::uint64_t>(p)));
        const bool itm = (type == OptionType::Call) ? (ST > K) : (ST < K);
        if (itm) sum += (kind == DigitalType::CashOrNothing) ? cash : ST;
    }
    return disc * sum / static_cast<double>(n_paths);
}

}  // namespace pricer

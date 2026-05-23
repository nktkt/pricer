// pricer/implied_vol.hpp — recover the volatility implied by a market price.
//
// The Black–Scholes price is strictly increasing in volatility, so the inverse
// is well defined. This uses a safeguarded Newton iteration (fast, via vega)
// that falls back to bisection whenever a step would leave the bracket — robust
// even for deep in/out-of-the-money quotes. It is the basic building block of
// volatility-surface construction and model calibration.
#pragma once
#include <cmath>
#include "pricer/black_scholes.hpp"

namespace pricer {

// Volatility such that black_scholes_price(type, S, K, r, sigma, T, q) == price.
inline double implied_vol(OptionType type, double price, double S, double K, double r, double T,
                          double tol = 1e-8, int max_iter = 100, double q = 0.0) {
    double lo = 1e-8, hi = 5.0;     // volatility bracket (price is monotone in sigma)
    double sigma = 0.2;             // initial guess
    const double sqrtT = std::sqrt(T);
    const double discq = std::exp(-q * T);

    for (int i = 0; i < max_iter; ++i) {
        const double p = black_scholes_price(type, S, K, r, sigma, T, q);
        const double diff = p - price;
        if (diff > 0) hi = sigma; else lo = sigma;   // keep the root bracketed
        if (std::fabs(diff) < tol) return sigma;

        const double d1 = (std::log(S / K) + (r - q + 0.5 * sigma * sigma) * T) / (sigma * sqrtT);
        const double vega = S * discq * norm_pdf(d1) * sqrtT;

        double next = sigma - diff / vega;            // Newton step
        if (!(next > lo && next < hi) || vega < 1e-12)
            next = 0.5 * (lo + hi);                    // safeguard: bisection step
        sigma = next;
    }
    return sigma;
}

}  // namespace pricer

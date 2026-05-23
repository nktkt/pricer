// pricer/bachelier.hpp — the Bachelier (normal) model.
//
// Black–Scholes models the underlying as a geometric Brownian motion, so prices
// stay positive — wrong for interest rates and spreads, which can go negative.
// The Bachelier (1900) model instead uses *arithmetic* Brownian motion on the
// forward,
//     F_T = F + sigma_n * sqrt(T) * Z,   Z ~ N(0,1),
// where sigma_n is an absolute ("normal") volatility in price units. The terminal
// forward is normal, so it can be negative, and the option has a clean closed
// form. This is the standard quote convention for rates options in a low-/
// negative-rate world.
//
// Everything is expressed on the forward F with an explicit discount factor `df`
// (e.g. e^{-rT}); the risk-free rate never appears on its own.
#pragma once
#include <algorithm>
#include <cmath>
#include <cstdint>

#include "pricer/black_scholes.hpp"  // OptionType, Greeks, norm_cdf, norm_pdf
#include "pricer/rng.hpp"            // cb_normal

namespace pricer {

// Bachelier price of a European option on a forward. `sigma_n` is the absolute
// (normal) volatility; `df` is the discount factor to expiry. At the money the
// price is simply df * sigma_n * sqrt(T) * phi(0); as sigma_n -> 0 it tends to the
// discounted intrinsic df * max(F - K, 0).
inline double bachelier_price(OptionType type, double F, double K, double sigma_n, double T,
                              double df = 1.0) {
    const double v = sigma_n * std::sqrt(T);  // total normal stdev over [0, T]
    if (v < 1e-300)
        return df * (type == OptionType::Call ? std::max(F - K, 0.0) : std::max(K - F, 0.0));
    const double d = (F - K) / v;
    if (type == OptionType::Call)
        return df * ((F - K) * norm_cdf(d) + v * norm_pdf(d));
    return df * ((K - F) * norm_cdf(-d) + v * norm_pdf(d));
}

// Greeks under the Bachelier model (sensitivities to the forward F and the normal
// vol sigma_n). delta = dP/dF, gamma = d2P/dF2, vega = dP/dsigma_n, theta is the
// time decay at a fixed discount factor (-dP/dT). `rho` is left 0 — the rate
// enters only through the exogenous `df`, so it is not modelled here.
inline Greeks bachelier_greeks(OptionType type, double F, double K, double sigma_n, double T,
                               double df = 1.0) {
    const double sqrtT = std::sqrt(T);
    const double v = sigma_n * sqrtT;
    const double d = (F - K) / v;
    const double pdf = norm_pdf(d);
    Greeks g{};
    g.price = bachelier_price(type, F, K, sigma_n, T, df);
    g.delta = (type == OptionType::Call) ? df * norm_cdf(d) : -df * norm_cdf(-d);
    g.gamma = df * pdf / v;
    g.vega = df * sqrtT * pdf;
    g.theta = -df * sigma_n * pdf / (2.0 * sqrtT);
    g.rho = 0.0;
    return g;
}

// Normal (Bachelier) implied volatility from a price, by safeguarded Newton. The
// at-the-money inversion sigma_n = price / (df * sqrt(T) * phi(0)) is exact at
// F = K and a strong starting guess elsewhere; vega > 0 makes Newton converge.
inline double bachelier_implied_vol(OptionType type, double price, double F, double K, double T,
                                    double df = 1.0, double tol = 1e-10, int max_iter = 100) {
    const double sqrtT = std::sqrt(T);
    constexpr double phi0 = 0.3989422804014327;  // 1 / sqrt(2*pi)
    double sigma = (price / df) / (sqrtT * phi0);  // exact ATM, good guess otherwise
    if (!(sigma > 0.0)) sigma = 1e-4;
    for (int i = 0; i < max_iter; ++i) {
        const double model = bachelier_price(type, F, K, sigma, T, df);
        const double diff = model - price;
        if (std::fabs(diff) < tol) break;
        const double d = (F - K) / (sigma * sqrtT);
        const double vega = df * sqrtT * norm_pdf(d);
        if (vega < 1e-14) break;
        sigma -= diff / vega;
        if (sigma <= 0.0) sigma = 1e-8;
    }
    return sigma;
}

// Monte Carlo price under arithmetic Brownian motion, with the counter-based RNG
// (reproducible). The terminal forward can be negative, which is the whole point
// of the normal model.
inline double bachelier_price_mc(OptionType type, double F, double K, double sigma_n, double T,
                                 long n_paths, std::uint64_t seed = 12345, double df = 1.0) {
    if (n_paths < 1) n_paths = 1;
    const double v = sigma_n * std::sqrt(T);
    double sum = 0.0;
    for (long p = 0; p < n_paths; ++p) {
        const double FT = F + v * cb_normal(seed, static_cast<std::uint64_t>(p));
        sum += (type == OptionType::Call) ? std::max(FT - K, 0.0) : std::max(K - FT, 0.0);
    }
    return df * sum / static_cast<double>(n_paths);
}

}  // namespace pricer

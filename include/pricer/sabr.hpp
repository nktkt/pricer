// pricer/sabr.hpp — the SABR stochastic-volatility model.
//
// SABR (Hagan, Kumar, Lesniak, Woodward 2002) models a forward F and its
// stochastic volatility α as
//     dF = α F^β dW1,   dα = ν α dW2,   dW1·dW2 = ρ dt
// with parameters α (initial vol level), β (CEV elasticity, fixed per market),
// ρ (spot/vol correlation, the skew) and ν (vol-of-vol, the smile convexity).
// It is the industry-standard model for interest-rate and FX smiles because
// Hagan's singular-perturbation result gives the Black (lognormal) implied
// volatility as a *closed form* — so a whole smile is one cheap formula, and
// calibration fits only α, ρ, ν at a chosen β.
//
// This header provides that implied-vol formula, option pricing via Black-76,
// least-squares calibration (through the shared Levenberg–Marquardt solver), and
// a Monte Carlo engine that simulates the SABR SDE directly to cross-check the
// asymptotic formula.
#pragma once
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <vector>

#include "pricer/black_scholes.hpp"  // OptionType, norm_cdf (via normal.hpp)
#include "pricer/optimize.hpp"
#include "pricer/rng.hpp"  // cb_normal

namespace pricer {

struct SabrParams {
    double alpha;  // initial volatility level (> 0)
    double beta;   // CEV elasticity in [0, 1] (fixed per market)
    double rho;    // forward/vol correlation in (-1, 1) (the skew)
    double nu;     // volatility of volatility (>= 0) (the smile)
};

// Hagan 2002 lognormal (Black) implied volatility for strike K and forward F at
// expiry T. Reduces to a flat smile σ = α when β = 1 and ν = 0, and the general
// branch matches the at-the-money branch continuously as K → F.
inline double sabr_implied_vol(double F, double K, double T, const SabrParams& p) {
    const double alpha = p.alpha, beta = p.beta, nu = p.nu;
    const double rho = std::clamp(p.rho, -0.999999, 0.999999);
    const double beta1 = 1.0 - beta;

    // Time-correction factor (common to both branches; uses the relevant power of
    // the forward/strike). Computed per branch below with the right base.
    if (std::fabs(F - K) < 1e-12 * std::max(1.0, F)) {  // at-the-money branch
        const double fb = std::pow(F, beta1);  // F^(1-β)
        const double corr =
            1.0 + ((beta1 * beta1 / 24.0) * alpha * alpha / std::pow(F, 2.0 * beta1) +
                   0.25 * rho * beta * nu * alpha / fb + (2.0 - 3.0 * rho * rho) / 24.0 * nu * nu) *
                      T;
        return alpha / fb * corr;
    }

    const double logFK = std::log(F / K);
    const double FK = F * K;
    const double fkbeta = std::pow(FK, beta1 / 2.0);  // (FK)^((1-β)/2)
    const double z = (nu / alpha) * fkbeta * logFK;
    // x(z); ratio z/x(z) → 1 as z → 0.
    const double xz = std::log((std::sqrt(1.0 - 2.0 * rho * z + z * z) + z - rho) / (1.0 - rho));
    const double zx = (std::fabs(z) < 1e-12) ? 1.0 : z / xz;

    const double denom =
        fkbeta * (1.0 + (beta1 * beta1 / 24.0) * logFK * logFK +
                  (beta1 * beta1 * beta1 * beta1 / 1920.0) * logFK * logFK * logFK * logFK);
    const double corr =
        1.0 + ((beta1 * beta1 / 24.0) * alpha * alpha / std::pow(FK, beta1) +
               0.25 * rho * beta * nu * alpha / fkbeta + (2.0 - 3.0 * rho * rho) / 24.0 * nu * nu) *
                  T;
    return (alpha / denom) * zx * corr;
}

// Price a European option on the forward via Black-76 using the SABR implied vol.
// `df` is the discount factor to expiry (e.g. exp(-r·T)).
inline double sabr_black_price(OptionType type, double F, double K, double T, double df,
                               const SabrParams& p) {
    const double vol = sabr_implied_vol(F, K, T, p);
    const double s = vol * std::sqrt(T);
    const double d1 = (std::log(F / K) + 0.5 * s * s) / s;
    const double d2 = d1 - s;
    if (type == OptionType::Call) return df * (F * norm_cdf(d1) - K * norm_cdf(d2));
    return df * (K * norm_cdf(-d2) - F * norm_cdf(-d1));
}

struct SabrFit {
    SabrParams params;
    double rms_vol_error;  // root-mean-square implied-vol error across the quotes
};

// Calibrate (α, ρ, ν) at a fixed β to a smile of market implied vols by least
// squares on vol error (Levenberg–Marquardt). α and ν are kept positive (fit on
// magnitude) and ρ stays in (-1, 1) — the optimizer is free to wander and the
// canonical values are read back. Needs at least three quotes.
inline SabrFit calibrate_sabr(double F, double T, const std::vector<double>& strikes,
                              const std::vector<double>& market_vols, double beta,
                              SabrParams guess = {0.2, 0.5, 0.0, 0.3}) {
    if (strikes.size() != market_vols.size() || strikes.size() < 3)
        throw std::invalid_argument("calibrate_sabr: need >= 3 matching quotes");
    const std::size_t n = strikes.size();

    // Seed alpha from the at-the-money vol so the optimizer starts at the right
    // scale regardless of beta: to leading order sigma_ATM ≈ alpha / F^(1-beta),
    // so alpha ≈ sigma_ATM · F^(1-beta). (Using a fixed guess fails when beta < 1
    // makes alpha a CEV-level vol far from a lognormal one.)
    std::size_t atm = 0;
    for (std::size_t i = 1; i < n; ++i)
        if (std::fabs(strikes[i] - F) < std::fabs(strikes[atm] - F)) atm = i;
    guess.alpha = market_vols[atm] * std::pow(F, 1.0 - beta);

    auto unpack = [beta](const std::vector<double>& p) {
        return SabrParams{std::fabs(p[0]), beta, std::clamp(p[1], -0.999999, 0.999999), std::fabs(p[2])};
    };
    auto residual = [&](const std::vector<double>& p, std::vector<double>& out) {
        const SabrParams q = unpack(p);
        for (std::size_t i = 0; i < n; ++i)
            out[i] = sabr_implied_vol(F, strikes[i], T, q) - market_vols[i];
    };

    std::vector<double> p0 = {guess.alpha, guess.rho, guess.nu};
    const opt::LMResult res = opt::levenberg_marquardt(residual, p0, static_cast<int>(n));
    const SabrParams fit = unpack(res.params);
    return {fit, std::sqrt(res.cost / static_cast<double>(n))};
}

// Monte Carlo price of a European option on the forward by simulating the SABR
// SDE directly with the counter-based RNG (reproducible). The vol process α is
// lognormal, so it is stepped exactly; the forward uses a full-truncation Euler
// step (F floored at 0, which the β < 1 CEV diffusion can reach). Used to
// cross-check the asymptotic Hagan formula.
inline double sabr_price_mc(OptionType type, double F0, double K, double T, double df,
                            const SabrParams& p, int n_steps, long n_paths,
                            std::uint64_t seed = 12345) {
    if (n_steps < 1) n_steps = 1;
    if (n_paths < 1) n_paths = 1;
    const double dt = T / n_steps;
    const double sqdt = std::sqrt(dt);
    const double rho = std::clamp(p.rho, -0.999999, 0.999999);
    const double rho2 = std::sqrt(1.0 - rho * rho);

    double sum = 0.0;
    for (long path = 0; path < n_paths; ++path) {
        double F = F0, a = p.alpha;
        for (int t = 0; t < n_steps; ++t) {
            const std::uint64_t base = (static_cast<std::uint64_t>(path) * n_steps + t) * 2;
            const double z1 = cb_normal(seed, base);
            const double z2 = cb_normal(seed, base + 1);
            const double w2 = rho * z1 + rho2 * z2;  // correlated Brownian increment for α
            const double Fpos = std::max(F, 0.0);
            F = F + a * std::pow(Fpos, p.beta) * sqdt * z1;
            if (F < 0.0) F = 0.0;  // full truncation / absorption at zero
            a *= std::exp(p.nu * sqdt * w2 - 0.5 * p.nu * p.nu * dt);  // exact lognormal step
        }
        sum += (type == OptionType::Call) ? std::max(F - K, 0.0) : std::max(K - F, 0.0);
    }
    return df * sum / static_cast<double>(n_paths);
}

}  // namespace pricer

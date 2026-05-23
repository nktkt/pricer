// pricer/exotics.hpp — path-dependent exotic options (Asian, barrier, lookback).
//
// Unlike a vanilla European option, whose value depends only on the terminal
// spot, these products depend on the whole price path:
//   * an *Asian* option pays on the average spot over a set of monitoring dates;
//   * a *barrier* option is knocked in or out if the spot ever crosses a level;
//   * a *lookback* option pays on the path's running minimum or maximum.
//
// Each family is priced two ways that cross-check each other: a closed form that
// is exact (Asian, under geometric averaging) or exact under continuous
// monitoring (barrier, lookback), and a Monte Carlo engine that steps the path
// with the counter-based RNG (so results are reproducible regardless of how the
// work is split). For barriers and lookbacks the MC applies the
// Broadie–Glasserman–Kou continuity correction so the discretely-stepped
// simulation converges to the continuous-monitoring closed form.
//
// All prices carry a continuous dividend yield q (cost of carry b = r - q),
// defaulting to 0.
#pragma once
#include <algorithm>
#include <cmath>
#include <cstdint>

#include "pricer/black_scholes.hpp"  // OptionType, Greeks, black_scholes_price, norm_cdf
#include "pricer/rng.hpp"            // cb_normal

namespace pricer {

namespace detail {
// Vanilla intrinsic at spot x (kept local so this header needn't pull in american.hpp).
inline double exotic_intrinsic(OptionType type, double x, double K) {
    return type == OptionType::Call ? std::max(x - K, 0.0) : std::max(K - x, 0.0);
}
// Broadie–Glasserman–Kou continuity-correction constant, -zeta(1/2)/sqrt(2*pi).
inline constexpr double kBGK = 0.5826;
}  // namespace detail

// =====================================================================
// Asian options (average price, fixed strike)
// =====================================================================

enum class AverageType { Arithmetic, Geometric };

// Closed-form price of a *discretely-monitored geometric-average* Asian option,
// monitored at t_i = i*T/n for i = 1..n. The geometric average of GBM is itself
// lognormal, so ln G ~ N(M, V) and the price is an exact Black–Scholes-style
// formula. This both prices the geometric Asian and is the validation backbone
// for the Monte Carlo engine below. (n = 1 reduces to Black–Scholes.)
inline double geometric_asian_price(OptionType type, double S, double K, double r, double sigma,
                                    double T, int n, double q = 0.0) {
    if (n < 1) n = 1;
    const double dt = T / n;
    const double mu_t = dt * (n + 1) / 2.0;  // (1/n) Σ t_i, the mean monitoring time
    const double V = sigma * sigma * dt * (n + 1) * (2.0 * n + 1) / (6.0 * n);  // Var[ln G]
    const double M = std::log(S) + (r - q - 0.5 * sigma * sigma) * mu_t;        // E[ln G]
    const double sd = std::sqrt(V);
    const double EG = std::exp(M + 0.5 * V);  // E[G], the average's forward
    const double d1 = (M - std::log(K) + V) / sd;
    const double d2 = d1 - sd;
    const double disc = std::exp(-r * T);
    if (type == OptionType::Call) return disc * (EG * norm_cdf(d1) - K * norm_cdf(d2));
    return disc * (K * norm_cdf(-d2) - EG * norm_cdf(-d1));
}

// Monte Carlo price of an average-price (fixed-strike) Asian option, averaging
// the spot over n_steps monitoring dates t_i = i*T/n (i = 1..n_steps) with the
// counter-based RNG. Arithmetic averaging has no closed form (the sum of
// lognormals is not lognormal); geometric averaging cross-checks against
// geometric_asian_price() to floating-point + MC tolerance.
inline double asian_price_mc(OptionType type, AverageType avg, double S, double K, double r,
                             double sigma, double T, int n_steps, long n_paths,
                             std::uint64_t seed = 12345, double q = 0.0) {
    if (n_steps < 1) n_steps = 1;
    if (n_paths < 1) n_paths = 1;
    const double dt = T / n_steps;
    const double drift = (r - q - 0.5 * sigma * sigma) * dt;
    const double vol = sigma * std::sqrt(dt);
    const double disc = std::exp(-r * T);

    double sum = 0.0;
    for (long p = 0; p < n_paths; ++p) {
        double s = S;
        double acc = 0.0;  // running sum (arithmetic) or sum of logs (geometric)
        for (int t = 1; t <= n_steps; ++t) {
            const std::uint64_t ctr = static_cast<std::uint64_t>(p) * n_steps + (t - 1);
            s *= std::exp(drift + vol * cb_normal(seed, ctr));
            acc += (avg == AverageType::Arithmetic) ? s : std::log(s);
        }
        const double avg_price =
            (avg == AverageType::Arithmetic) ? acc / n_steps : std::exp(acc / n_steps);
        sum += detail::exotic_intrinsic(type, avg_price, K);
    }
    return disc * sum / static_cast<double>(n_paths);
}

// =====================================================================
// Barrier options (single barrier, continuous monitoring, no rebate)
// =====================================================================

enum class BarrierType { UpOut, UpIn, DownOut, DownIn };

namespace detail {
// The four Reiner–Rubinstein building blocks A, B, C, D for standard barrier
// options, parameterized by phi (call=+1, put=-1) and eta (down=+1, up=-1).
struct RRBlocks {
    double A, B, C, D;
};

inline RRBlocks rr_blocks(double S, double K, double Bar, double r, double sigma, double T,
                          double q, double phi, double eta) {
    const double b = r - q;  // cost of carry
    const double vsqrtT = sigma * std::sqrt(T);
    const double mu = (b - 0.5 * sigma * sigma) / (sigma * sigma);
    const double x1 = std::log(S / K) / vsqrtT + (1.0 + mu) * vsqrtT;
    const double x2 = std::log(S / Bar) / vsqrtT + (1.0 + mu) * vsqrtT;
    const double y1 = std::log(Bar * Bar / (S * K)) / vsqrtT + (1.0 + mu) * vsqrtT;
    const double y2 = std::log(Bar / S) / vsqrtT + (1.0 + mu) * vsqrtT;
    const double ebrT = std::exp((b - r) * T);
    const double erT = std::exp(-r * T);
    const double pm1 = std::pow(Bar / S, 2.0 * (mu + 1.0));
    const double pm = std::pow(Bar / S, 2.0 * mu);
    RRBlocks o;
    o.A = phi * S * ebrT * norm_cdf(phi * x1) - phi * K * erT * norm_cdf(phi * x1 - phi * vsqrtT);
    o.B = phi * S * ebrT * norm_cdf(phi * x2) - phi * K * erT * norm_cdf(phi * x2 - phi * vsqrtT);
    o.C = phi * S * ebrT * pm1 * norm_cdf(eta * y1) -
          phi * K * erT * pm * norm_cdf(eta * y1 - eta * vsqrtT);
    o.D = phi * S * ebrT * pm1 * norm_cdf(eta * y2) -
          phi * K * erT * pm * norm_cdf(eta * y2 - eta * vsqrtT);
    return o;
}
}  // namespace detail

// Closed-form price of a standard single-barrier option (continuous monitoring,
// no rebate) via the Reiner–Rubinstein formulas. By construction a knock-in plus
// the matching knock-out equals the vanilla option (parity): exactly one is alive
// at expiry. Inputs are assumed not already knocked at t=0; if they are, an
// already-out option is worth 0 and an already-in option is the vanilla.
inline double barrier_price(OptionType type, BarrierType barrier, double S, double K, double Bar,
                            double r, double sigma, double T, double q = 0.0) {
    const bool up = (barrier == BarrierType::UpOut || barrier == BarrierType::UpIn);
    const bool is_out = (barrier == BarrierType::UpOut || barrier == BarrierType::DownOut);
    const bool already = up ? (S >= Bar) : (S <= Bar);
    if (already) return is_out ? 0.0 : black_scholes_price(type, S, K, r, sigma, T, q);

    const double phi = (type == OptionType::Call) ? 1.0 : -1.0;
    const double eta = up ? -1.0 : 1.0;
    const detail::RRBlocks k = detail::rr_blocks(S, K, Bar, r, sigma, T, q, phi, eta);
    const bool KgeB = (K >= Bar);
    switch (barrier) {
        case BarrierType::DownIn:
            return type == OptionType::Call ? (KgeB ? k.C : k.A - k.B + k.D)
                                            : (KgeB ? k.B - k.C + k.D : k.A);
        case BarrierType::UpIn:
            return type == OptionType::Call ? (KgeB ? k.A : k.B - k.C + k.D)
                                            : (KgeB ? k.A - k.B + k.D : k.C);
        case BarrierType::DownOut:
            return type == OptionType::Call ? (KgeB ? k.A - k.C : k.B - k.D)
                                            : (KgeB ? k.A - k.B + k.C - k.D : 0.0);
        case BarrierType::UpOut:
            return type == OptionType::Call ? (KgeB ? 0.0 : k.A - k.B + k.C - k.D)
                                            : (KgeB ? k.B - k.D : k.A - k.C);
    }
    return 0.0;
}

// Monte Carlo price of a single-barrier option, stepping the path over n_steps
// monitoring dates with the counter-based RNG. By default applies the
// Broadie–Glasserman–Kou continuity correction — shifting the monitored barrier
// inward by exp(∓0.5826 σ√dt) (up barriers down, down barriers up) so the
// discrete grid catches the crossings continuous monitoring would — making the
// discretely-stepped path converge to the continuous-monitoring price returned
// by barrier_price(). Set continuity_correction=false to price a genuinely
// discretely-monitored barrier.
inline double barrier_price_mc(OptionType type, BarrierType barrier, double S, double K, double Bar,
                               double r, double sigma, double T, int n_steps, long n_paths,
                               std::uint64_t seed = 12345, double q = 0.0,
                               bool continuity_correction = true) {
    if (n_steps < 1) n_steps = 1;
    if (n_paths < 1) n_paths = 1;
    const double dt = T / n_steps;
    const double drift = (r - q - 0.5 * sigma * sigma) * dt;
    const double vol = sigma * std::sqrt(dt);
    const double disc = std::exp(-r * T);
    const bool up = (barrier == BarrierType::UpOut || barrier == BarrierType::UpIn);
    const bool is_out = (barrier == BarrierType::UpOut || barrier == BarrierType::DownOut);

    double Bm = Bar;
    if (continuity_correction)
        Bm = Bar * std::exp((up ? -1.0 : 1.0) * detail::kBGK * sigma * std::sqrt(dt));

    double sum = 0.0;
    for (long p = 0; p < n_paths; ++p) {
        double s = S;
        bool hit = false;
        for (int t = 1; t <= n_steps; ++t) {
            const std::uint64_t ctr = static_cast<std::uint64_t>(p) * n_steps + (t - 1);
            s *= std::exp(drift + vol * cb_normal(seed, ctr));
            if (up ? (s >= Bm) : (s <= Bm)) hit = true;
        }
        const bool alive = is_out ? !hit : hit;
        if (alive) sum += detail::exotic_intrinsic(type, s, K);
    }
    return disc * sum / static_cast<double>(n_paths);
}

// =====================================================================
// Lookback options (floating strike, continuous monitoring, freshly issued)
// =====================================================================

// Closed-form floating-strike lookback (Conze–Viswanathan / Goldman–Sosin–Gatto)
// for a freshly issued option whose running extreme equals the spot. The call
// pays S_T − min S (buy at the lowest price seen), the put pays max S − S_T (sell
// at the highest). Cost of carry b = r − q must be non-zero (the formula has a
// 1/(2b) term).
inline double lookback_floating_price(OptionType type, double S, double r, double sigma, double T,
                                      double q = 0.0) {
    const double b = r - q;
    const double sqrtT = std::sqrt(T);
    const double vsqrtT = sigma * sqrtT;
    const double ebrT = std::exp((b - r) * T);
    const double erT = std::exp(-r * T);
    const double kk = sigma * sigma / (2.0 * b);
    if (type == OptionType::Call) {
        const double a1 = (b + 0.5 * sigma * sigma) * T / vsqrtT;  // ln(S/min)=0 (fresh)
        const double a2 = a1 - vsqrtT;
        return S * ebrT * norm_cdf(a1) - S * erT * norm_cdf(a2) +
               S * erT * kk * (norm_cdf(-a1 + (2.0 * b / sigma) * sqrtT) -
                               std::exp(b * T) * norm_cdf(-a1));
    }
    const double b1 = (b + 0.5 * sigma * sigma) * T / vsqrtT;  // ln(S/max)=0 (fresh)
    const double b2 = b1 - vsqrtT;
    return S * erT * norm_cdf(-b2) - S * ebrT * norm_cdf(-b1) +
           S * erT * kk * (-norm_cdf(b1 - (2.0 * b / sigma) * sqrtT) +
                           std::exp(b * T) * norm_cdf(b1));
}

// Monte Carlo price of a freshly-issued floating-strike lookback, tracking the
// path's running min/max over n_steps dates with the counter-based RNG. The
// discrete max/min understates the continuous extreme, so by default the
// Broadie–Glasserman–Kou correction nudges the running min down / max up by
// exp(±0.5826 σ√dt) to converge to lookback_floating_price().
inline double lookback_floating_price_mc(OptionType type, double S, double r, double sigma, double T,
                                         int n_steps, long n_paths, std::uint64_t seed = 12345,
                                         double q = 0.0, bool continuity_correction = true) {
    if (n_steps < 1) n_steps = 1;
    if (n_paths < 1) n_paths = 1;
    const double dt = T / n_steps;
    const double drift = (r - q - 0.5 * sigma * sigma) * dt;
    const double vol = sigma * std::sqrt(dt);
    const double disc = std::exp(-r * T);
    const double corr = continuity_correction ? std::exp(detail::kBGK * sigma * std::sqrt(dt)) : 1.0;

    double sum = 0.0;
    for (long p = 0; p < n_paths; ++p) {
        double s = S, mn = S, mx = S;
        for (int t = 1; t <= n_steps; ++t) {
            const std::uint64_t ctr = static_cast<std::uint64_t>(p) * n_steps + (t - 1);
            s *= std::exp(drift + vol * cb_normal(seed, ctr));
            mn = std::min(mn, s);
            mx = std::max(mx, s);
        }
        mn /= corr;
        mx *= corr;
        sum += (type == OptionType::Call) ? (s - mn) : (mx - s);
    }
    return disc * sum / static_cast<double>(n_paths);
}

}  // namespace pricer

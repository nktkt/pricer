// pricer/basket.hpp — multi-asset options on correlated geometric Brownian motion.
//
// These products pay on several underlyings at once, so pricing them means
// simulating correlated GBM: draw independent normals, correlate them through the
// Cholesky factor of the correlation matrix, and evolve each asset. The header
// covers the two classic families, each with a closed form that anchors a Monte
// Carlo engine:
//   * basket options on a weighted portfolio Σ w_i S_i — the geometric-average
//     basket is lognormal (exact closed form); the arithmetic basket (the real
//     product) has no closed form and is priced by MC, with the geometric basket
//     as a cross-check;
//   * spread options max(S1 − S2 − K, 0) — at K = 0 this is an exchange option
//     with Margrabe's exact formula; for K ≠ 0 there is Kirk's approximation,
//     cross-checked against MC.
//
// Terminal payoffs need only the assets at expiry, so each path is a single exact
// GBM jump (no time grid). The counter-based RNG keeps results reproducible.
#pragma once
#include <cmath>
#include <cstdint>
#include <vector>

#include "pricer/black_scholes.hpp"  // OptionType, norm_cdf
#include "pricer/exotics.hpp"        // AverageType, detail::exotic_intrinsic
#include "pricer/rng.hpp"            // cb_normal

namespace pricer {

namespace detail {
// Lower-triangular Cholesky factor L of a symmetric positive-(semi)definite
// matrix A, so that L·Lᵀ = A. Diagonal negatives are clamped to 0, which keeps a
// borderline correlation matrix usable.
inline std::vector<std::vector<double>> cholesky(const std::vector<std::vector<double>>& A) {
    const std::size_t n = A.size();
    std::vector<std::vector<double>> L(n, std::vector<double>(n, 0.0));
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j <= i; ++j) {
            double s = A[i][j];
            for (std::size_t k = 0; k < j; ++k) s -= L[i][k] * L[j][k];
            if (i == j)
                L[i][j] = std::sqrt(s > 0.0 ? s : 0.0);
            else
                L[i][j] = (L[j][j] > 0.0) ? s / L[j][j] : 0.0;
        }
    return L;
}

// Dividend yield of asset i (0 if the q vector is omitted).
inline double q_of(const std::vector<double>& q, std::size_t i) { return q.empty() ? 0.0 : q[i]; }
}  // namespace detail

// =====================================================================
// Basket options (weighted portfolio of assets)
// =====================================================================

// Closed-form price of a *geometric-average* basket option, payoff on
// G = Π_i S_i(T)^{w_i}. ln G is normal, so this is an exact Black–Scholes-style
// formula and the validation backbone for the Monte Carlo engine. With one asset
// and weight 1 it reduces to Black–Scholes.
inline double geometric_basket_price(OptionType type, const std::vector<double>& S,
                                     const std::vector<double>& w, double K, double r,
                                     const std::vector<double>& sigma,
                                     const std::vector<std::vector<double>>& corr, double T,
                                     const std::vector<double>& q = {}) {
    const std::size_t n = S.size();
    double M = 0.0, V = 0.0;  // mean and variance of ln G
    for (std::size_t i = 0; i < n; ++i)
        M += w[i] * (std::log(S[i]) + (r - detail::q_of(q, i) - 0.5 * sigma[i] * sigma[i]) * T);
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < n; ++j)
            V += w[i] * w[j] * sigma[i] * sigma[j] * corr[i][j] * T;
    const double sd = std::sqrt(V);
    const double EG = std::exp(M + 0.5 * V);
    const double d1 = (M - std::log(K) + V) / sd;
    const double d2 = d1 - sd;
    const double disc = std::exp(-r * T);
    if (type == OptionType::Call) return disc * (EG * norm_cdf(d1) - K * norm_cdf(d2));
    return disc * (K * norm_cdf(-d2) - EG * norm_cdf(-d1));
}

// Monte Carlo price of a basket option (arithmetic Σ w_i S_i(T), or geometric
// Π S_i(T)^{w_i}) under correlated GBM. Correlated normals come from the Cholesky
// factor of `corr`. The arithmetic basket has no closed form; the geometric mode
// cross-checks against geometric_basket_price().
inline double basket_price_mc(OptionType type, AverageType avg, const std::vector<double>& S,
                              const std::vector<double>& w, double K, double r,
                              const std::vector<double>& sigma,
                              const std::vector<std::vector<double>>& corr, double T, long n_paths,
                              std::uint64_t seed = 12345, const std::vector<double>& q = {}) {
    if (n_paths < 1) n_paths = 1;
    const std::size_t n = S.size();
    const auto L = detail::cholesky(corr);
    std::vector<double> drift(n), vol(n);
    for (std::size_t i = 0; i < n; ++i) {
        drift[i] = (r - detail::q_of(q, i) - 0.5 * sigma[i] * sigma[i]) * T;
        vol[i] = sigma[i] * std::sqrt(T);
    }
    const double disc = std::exp(-r * T);

    double sum = 0.0;
    std::vector<double> Z(n);
    for (long p = 0; p < n_paths; ++p) {
        for (std::size_t k = 0; k < n; ++k)
            Z[k] = cb_normal(seed, static_cast<std::uint64_t>(p) * n + k);
        double basket = 0.0, lg = 0.0;
        for (std::size_t i = 0; i < n; ++i) {
            double x = 0.0;  // correlated normal X_i = (L·Z)_i
            for (std::size_t k = 0; k <= i; ++k) x += L[i][k] * Z[k];
            const double Si = S[i] * std::exp(drift[i] + vol[i] * x);
            if (avg == AverageType::Arithmetic)
                basket += w[i] * Si;
            else
                lg += w[i] * std::log(Si);
        }
        if (avg == AverageType::Geometric) basket = std::exp(lg);
        sum += detail::exotic_intrinsic(type, basket, K);
    }
    return disc * sum / static_cast<double>(n_paths);
}

// =====================================================================
// Spread / exchange options (two assets)
// =====================================================================

// Margrabe's exact price of a European exchange option: the right to give up
// asset 2 and receive asset 1 at expiry, payoff max(S1(T) − S2(T), 0). The
// risk-free rate cancels (the strike is itself an asset), leaving only the two
// dividend yields and the spread volatility √(σ1² + σ2² − 2ρσ1σ2).
inline double margrabe_exchange_price(double S1, double S2, double sigma1, double sigma2, double rho,
                                      double T, double q1 = 0.0, double q2 = 0.0) {
    const double s = std::sqrt(sigma1 * sigma1 + sigma2 * sigma2 - 2.0 * rho * sigma1 * sigma2);
    const double sT = s * std::sqrt(T);
    const double d1 = (std::log(S1 / S2) + (q2 - q1 + 0.5 * s * s) * T) / sT;
    const double d2 = d1 - sT;
    return S1 * std::exp(-q1 * T) * norm_cdf(d1) - S2 * std::exp(-q2 * T) * norm_cdf(d2);
}

// Kirk's approximation for a spread option max(S1(T) − S2(T) − K, 0): treats
// S2 + K as a single lognormal and applies a Margrabe-style formula on the
// forwards. Exact at K = 0 (it reduces to Margrabe); for K ≠ 0 it is a very good
// approximation, cross-checked against spread_price_mc().
inline double spread_kirk_price(OptionType type, double S1, double S2, double K, double r,
                                double sigma1, double sigma2, double rho, double T, double q1 = 0.0,
                                double q2 = 0.0) {
    const double F1 = S1 * std::exp((r - q1) * T);
    const double F2 = S2 * std::exp((r - q2) * T);
    const double a = F2 / (F2 + K);
    const double s = std::sqrt(sigma1 * sigma1 - 2.0 * rho * sigma1 * sigma2 * a + sigma2 * sigma2 * a * a);
    const double sT = s * std::sqrt(T);
    const double d1 = (std::log(F1 / (F2 + K)) + 0.5 * s * s * T) / sT;
    const double d2 = d1 - sT;
    const double disc = std::exp(-r * T);
    if (type == OptionType::Call) return disc * (F1 * norm_cdf(d1) - (F2 + K) * norm_cdf(d2));
    return disc * ((F2 + K) * norm_cdf(-d2) - F1 * norm_cdf(-d1));
}

// Monte Carlo price of a spread option max(S1(T) − S2(T) − K, 0) (call) under
// correlated two-asset GBM. Matches Margrabe exactly at K = 0 and Kirk closely
// otherwise.
inline double spread_price_mc(OptionType type, double S1, double S2, double K, double r,
                              double sigma1, double sigma2, double rho, double T, long n_paths,
                              std::uint64_t seed = 12345, double q1 = 0.0, double q2 = 0.0) {
    if (n_paths < 1) n_paths = 1;
    const double dr1 = (r - q1 - 0.5 * sigma1 * sigma1) * T, v1 = sigma1 * std::sqrt(T);
    const double dr2 = (r - q2 - 0.5 * sigma2 * sigma2) * T, v2 = sigma2 * std::sqrt(T);
    const double rho2 = std::sqrt(1.0 - rho * rho);
    const double disc = std::exp(-r * T);
    double sum = 0.0;
    for (long p = 0; p < n_paths; ++p) {
        const double z1 = cb_normal(seed, static_cast<std::uint64_t>(p) * 2);
        const double z2 = cb_normal(seed, static_cast<std::uint64_t>(p) * 2 + 1);
        const double ST1 = S1 * std::exp(dr1 + v1 * z1);
        const double ST2 = S2 * std::exp(dr2 + v2 * (rho * z1 + rho2 * z2));
        sum += detail::exotic_intrinsic(type, ST1 - ST2, K);
    }
    return disc * sum / static_cast<double>(n_paths);
}

}  // namespace pricer

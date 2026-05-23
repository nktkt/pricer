// pricer/rainbow.hpp — rainbow (two-asset best-of / worst-of) options.
//
// A rainbow option pays on the maximum ("best-of") or minimum ("worst-of") of two
// correlated assets, e.g. max(max(S1,S2) − K, 0). Stulz (1982) gives exact
// closed forms for these in terms of the *bivariate* normal CDF, so this header
// first implements that CDF (Drezner's 1978 reduction to a 2-D Gauss quadrature)
// and then the option formulas. The puts come from the calls by put–call parity,
// using the discounted forward on the max/min (which is itself a vanilla position
// plus a Margrabe exchange option). A correlated two-asset Monte Carlo engine
// cross-checks everything.
//
// The cleanest validation is a parity that pins both call formulas at once: a call
// on the max plus a call on the min equals two vanilla calls,
//   C_max(K) + C_min(K) = c(S1,K) + c(S2,K).
#pragma once
#include <algorithm>
#include <cmath>
#include <cstdint>

#include "pricer/basket.hpp"         // margrabe_exchange_price, detail::exotic_intrinsic
#include "pricer/black_scholes.hpp"  // OptionType, norm_cdf, black_scholes_call
#include "pricer/rng.hpp"            // cb_normal

namespace pricer {

// Standard bivariate normal CDF M(a, b; rho) = P(X <= a, Y <= b) for a standard
// bivariate normal with correlation rho. Drezner (1978): a 4×4 Gauss quadrature
// on the negative-orthant case, with sign reflections reducing the others to it.
inline double bivariate_normal_cdf(double a, double b, double rho) {
    constexpr double kPi = 3.14159265358979323846;
    auto sgn = [](double x) { return x >= 0.0 ? 1.0 : -1.0; };

    if (a <= 0.0 && b <= 0.0 && rho <= 0.0) {
        const double denom = std::sqrt(2.0 * (1.0 - rho * rho));
        const double ap = a / denom, bp = b / denom;
        static const double X[4] = {0.3253030, 0.4211071, 0.1334425, 0.006374323};
        static const double Y[4] = {0.1337764, 0.6243247, 1.3425378, 2.2626645};
        double sum = 0.0;
        for (int i = 0; i < 4; ++i)
            for (int j = 0; j < 4; ++j)
                sum += X[i] * X[j] * std::exp(ap * (2.0 * Y[i] - ap) + bp * (2.0 * Y[j] - bp) +
                                              2.0 * rho * (Y[i] - ap) * (Y[j] - bp));
        return std::sqrt(1.0 - rho * rho) / kPi * sum;
    }
    if (a * b * rho <= 0.0) {
        if (a <= 0.0 && b >= 0.0 && rho >= 0.0) return norm_cdf(a) - bivariate_normal_cdf(a, -b, -rho);
        if (a >= 0.0 && b <= 0.0 && rho >= 0.0) return norm_cdf(b) - bivariate_normal_cdf(-a, b, -rho);
        if (a >= 0.0 && b >= 0.0 && rho <= 0.0)
            return norm_cdf(a) + norm_cdf(b) - 1.0 + bivariate_normal_cdf(-a, -b, rho);
    }
    const double denum = std::sqrt(a * a - 2.0 * rho * a * b + b * b);
    const double rho1 = (rho * a - b) * sgn(a) / denum;
    const double rho2 = (rho * b - a) * sgn(b) / denum;
    const double delta = (1.0 - sgn(a) * sgn(b)) / 4.0;
    return bivariate_normal_cdf(a, 0.0, rho1) + bivariate_normal_cdf(b, 0.0, rho2) - delta;
}

enum class RainbowType { Max, Min };  // option on the best-of (max) or worst-of (min)

namespace detail {
// Stulz call on the maximum of two assets, payoff max(max(S1,S2) − K, 0).
inline double stulz_call_max(double S1, double S2, double K, double r, double sigma1, double sigma2,
                             double rho, double T, double q1, double q2) {
    const double b1 = r - q1, b2 = r - q2;
    const double sig = std::sqrt(sigma1 * sigma1 + sigma2 * sigma2 - 2.0 * rho * sigma1 * sigma2);
    const double sT = std::sqrt(T);
    const double d = (std::log(S1 / S2) + (b1 - b2 + 0.5 * sig * sig) * T) / (sig * sT);
    const double y1 = (std::log(S1 / K) + (b1 + 0.5 * sigma1 * sigma1) * T) / (sigma1 * sT);
    const double y2 = (std::log(S2 / K) + (b2 + 0.5 * sigma2 * sigma2) * T) / (sigma2 * sT);
    const double rho1 = (sigma1 - rho * sigma2) / sig;
    const double rho2 = (sigma2 - rho * sigma1) / sig;
    return S1 * std::exp((b1 - r) * T) * bivariate_normal_cdf(y1, d, rho1) +
           S2 * std::exp((b2 - r) * T) * bivariate_normal_cdf(y2, -d + sig * sT, rho2) -
           K * std::exp(-r * T) *
               (1.0 - bivariate_normal_cdf(-y1 + sigma1 * sT, -y2 + sigma2 * sT, rho));
}

// Stulz call on the minimum of two assets, payoff max(min(S1,S2) − K, 0).
inline double stulz_call_min(double S1, double S2, double K, double r, double sigma1, double sigma2,
                             double rho, double T, double q1, double q2) {
    const double b1 = r - q1, b2 = r - q2;
    const double sig = std::sqrt(sigma1 * sigma1 + sigma2 * sigma2 - 2.0 * rho * sigma1 * sigma2);
    const double sT = std::sqrt(T);
    const double d = (std::log(S1 / S2) + (b1 - b2 + 0.5 * sig * sig) * T) / (sig * sT);
    const double y1 = (std::log(S1 / K) + (b1 + 0.5 * sigma1 * sigma1) * T) / (sigma1 * sT);
    const double y2 = (std::log(S2 / K) + (b2 + 0.5 * sigma2 * sigma2) * T) / (sigma2 * sT);
    const double rho1 = (sigma1 - rho * sigma2) / sig;
    const double rho2 = (sigma2 - rho * sigma1) / sig;
    return S1 * std::exp((b1 - r) * T) * bivariate_normal_cdf(y1, -d, -rho1) +
           S2 * std::exp((b2 - r) * T) * bivariate_normal_cdf(y2, d - sig * sT, -rho2) -
           K * std::exp(-r * T) * bivariate_normal_cdf(y1 - sigma1 * sT, y2 - sigma2 * sT, rho);
}
}  // namespace detail

// Closed-form price of a rainbow option on two assets (Stulz). The call formulas
// are exact; the puts come from put–call parity, P = C − F + K·e^{-rT}, where F is
// the discounted forward on the max/min: max(S1,S2) = S1 + max(S2−S1, 0) and
// min(S1,S2) = S1 − max(S1−S2, 0), so F is a vanilla leg plus/minus a Margrabe
// exchange option.
inline double rainbow_price(OptionType type, RainbowType which, double S1, double S2, double K,
                            double r, double sigma1, double sigma2, double rho, double T,
                            double q1 = 0.0, double q2 = 0.0) {
    const double call = (which == RainbowType::Max)
                            ? detail::stulz_call_max(S1, S2, K, r, sigma1, sigma2, rho, T, q1, q2)
                            : detail::stulz_call_min(S1, S2, K, r, sigma1, sigma2, rho, T, q1, q2);
    if (type == OptionType::Call) return call;
    // Put via parity using the discounted forward on the max / min.
    const double disc = std::exp(-r * T);
    const double fwd = (which == RainbowType::Max)
                           ? S1 * std::exp(-q1 * T) +
                                 margrabe_exchange_price(S2, S1, sigma2, sigma1, rho, T, q2, q1)
                           : S1 * std::exp(-q1 * T) -
                                 margrabe_exchange_price(S1, S2, sigma1, sigma2, rho, T, q1, q2);
    return call - fwd + K * disc;
}

// Monte Carlo price of a rainbow option under correlated two-asset GBM, with the
// counter-based RNG (reproducible). Handles call/put on the max/min directly.
inline double rainbow_price_mc(OptionType type, RainbowType which, double S1, double S2, double K,
                               double r, double sigma1, double sigma2, double rho, double T,
                               long n_paths, std::uint64_t seed = 12345, double q1 = 0.0,
                               double q2 = 0.0) {
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
        const double u = (which == RainbowType::Max) ? std::max(ST1, ST2) : std::min(ST1, ST2);
        sum += detail::exotic_intrinsic(type, u, K);
    }
    return disc * sum / static_cast<double>(n_paths);
}

}  // namespace pricer

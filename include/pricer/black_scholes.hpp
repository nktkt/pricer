// pricer/black_scholes.hpp — closed-form European option pricing and Greeks.
#pragma once
#include <cmath>
#include "pricer/normal.hpp"

namespace pricer {

enum class OptionType { Call, Put };

// Risk sensitivities of an option, returned together to avoid recomputing d1/d2.
struct Greeks {
    double price;  // option value
    double delta;  // d(price)/d(S)
    double gamma;  // d^2(price)/d(S)^2
    double vega;   // d(price)/d(sigma)  (per 1.0 of vol; divide by 100 for 1%)
    double theta;  // d(price)/d(t)      (per year)
    double rho;    // d(price)/d(r)      (per 1.0 of rate)
};

// Black–Scholes price of a European option (Merton form with a continuous
// dividend yield q; q = 0 is the classic no-dividend Black–Scholes).
//   S: spot, K: strike, r: risk-free rate, sigma: volatility, T: time to expiry
//   (years), q: continuous dividend yield
inline double black_scholes_price(OptionType type, double S, double K,
                                  double r, double sigma, double T, double q = 0.0) {
    const double sqrtT = std::sqrt(T);
    const double d1 = (std::log(S / K) + (r - q + 0.5 * sigma * sigma) * T) / (sigma * sqrtT);
    const double d2 = d1 - sigma * sqrtT;
    const double disc = std::exp(-r * T);
    const double discq = std::exp(-q * T);  // the spot is carried at r - q
    if (type == OptionType::Call)
        return S * discq * norm_cdf(d1) - K * disc * norm_cdf(d2);
    else
        return K * disc * norm_cdf(-d2) - S * discq * norm_cdf(-d1);
}

inline double black_scholes_call(double S, double K, double r, double sigma, double T,
                                 double q = 0.0) {
    return black_scholes_price(OptionType::Call, S, K, r, sigma, T, q);
}

inline double black_scholes_put(double S, double K, double r, double sigma, double T,
                                double q = 0.0) {
    return black_scholes_price(OptionType::Put, S, K, r, sigma, T, q);
}

// Closed-form Greeks for a European option (with continuous dividend yield q).
inline Greeks black_scholes_greeks(OptionType type, double S, double K,
                                   double r, double sigma, double T, double q = 0.0) {
    const double sqrtT = std::sqrt(T);
    const double d1 = (std::log(S / K) + (r - q + 0.5 * sigma * sigma) * T) / (sigma * sqrtT);
    const double d2 = d1 - sigma * sqrtT;
    const double disc = std::exp(-r * T);
    const double discq = std::exp(-q * T);
    const double pdf_d1 = norm_pdf(d1);

    Greeks g{};
    g.price = black_scholes_price(type, S, K, r, sigma, T, q);
    g.gamma = discq * pdf_d1 / (S * sigma * sqrtT);
    g.vega  = S * discq * pdf_d1 * sqrtT;
    if (type == OptionType::Call) {
        g.delta = discq * norm_cdf(d1);
        g.theta = -(S * discq * pdf_d1 * sigma) / (2 * sqrtT) - r * K * disc * norm_cdf(d2) +
                  q * S * discq * norm_cdf(d1);
        g.rho   = K * T * disc * norm_cdf(d2);
    } else {
        g.delta = discq * (norm_cdf(d1) - 1.0);
        g.theta = -(S * discq * pdf_d1 * sigma) / (2 * sqrtT) + r * K * disc * norm_cdf(-d2) -
                  q * S * discq * norm_cdf(-d1);
        g.rho   = -K * T * disc * norm_cdf(-d2);
    }
    return g;
}

}  // namespace pricer

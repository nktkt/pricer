// pricer/heston.hpp — the Heston stochastic-volatility model and its calibration.
//
// Prices European options under Heston via the semi-analytic characteristic-
// function formula (the numerically stable "little trap" form), integrated with
// Gauss–Legendre quadrature:
//     C = S·P1 − K·e^{-rT}·P2,   P_j = 1/2 + (1/π) ∫₀^∞ Re[e^{-iφ ln K} f_j(φ)/(iφ)] dφ
// `calibrate_heston` fits the five parameters to market option prices by least
// squares using the Levenberg–Marquardt solver.
#pragma once
#include <complex>
#include <cmath>
#include <vector>

#include "pricer/black_scholes.hpp"  // OptionType
#include "pricer/optimize.hpp"
#include "pricer/quadrature.hpp"

namespace pricer {

struct HestonParams {
    double kappa;  // mean-reversion speed
    double theta;  // long-run variance
    double sigma;  // vol of variance
    double rho;    // spot/variance correlation
    double v0;     // initial variance
};

namespace detail {
// Heston characteristic function term f_j (j = 1, 2), little-trap formulation.
inline std::complex<double> heston_cf(double phi, int j, const HestonParams& p, double S,
                                      double r, double T) {
    using C = std::complex<double>;
    const C i(0.0, 1.0);
    const double u = (j == 1) ? 0.5 : -0.5;
    const double b = (j == 1) ? (p.kappa - p.rho * p.sigma) : p.kappa;
    const C rsi = p.rho * p.sigma * i * phi;
    const C d = std::sqrt((rsi - b) * (rsi - b) - p.sigma * p.sigma * (2.0 * u * i * phi - phi * phi));
    const C g = (b - rsi - d) / (b - rsi + d);  // c = 1/g (trap form)
    const C edt = std::exp(-d * T);
    const C D = ((b - rsi - d) / (p.sigma * p.sigma)) * ((1.0 - edt) / (1.0 - g * edt));
    const C Cc = r * i * phi * T + (p.kappa * p.theta / (p.sigma * p.sigma)) *
                                       ((b - rsi - d) * T - 2.0 * std::log((1.0 - g * edt) / (1.0 - g)));
    return std::exp(Cc + D * p.v0 + i * phi * std::log(S));
}
}  // namespace detail

// Heston European call price.
inline double heston_call(const HestonParams& p, double S, double K, double r, double T) {
    using C = std::complex<double>;
    const C i(0.0, 1.0);
    static const auto gl = gauss_legendre(128, 1e-8, 200.0);  // fixed grid, computed once
    const auto& x = gl.first;
    const auto& w = gl.second;
    const double lnK = std::log(K);

    double P1 = 0.0, P2 = 0.0;
    for (std::size_t n = 0; n < x.size(); ++n) {
        const double phi = x[n];
        const C ekp = std::exp(-i * phi * lnK);
        P1 += w[n] * std::real(ekp * detail::heston_cf(phi, 1, p, S, r, T) / (i * phi));
        P2 += w[n] * std::real(ekp * detail::heston_cf(phi, 2, p, S, r, T) / (i * phi));
    }
    const double pi = 3.14159265358979323846;
    P1 = 0.5 + P1 / pi;
    P2 = 0.5 + P2 / pi;
    return S * P1 - K * std::exp(-r * T) * P2;
}

// Heston price for either option type (put via put–call parity).
inline double heston_price(OptionType type, const HestonParams& p, double S, double K, double r,
                           double T) {
    const double call = heston_call(p, S, K, r, T);
    return (type == OptionType::Call) ? call : call - S + K * std::exp(-r * T);
}

// Calibrate the five Heston parameters to market call prices at quotes
// (strikes[i], expiries[i]) by least squares (Levenberg–Marquardt).
inline HestonParams calibrate_heston(double S, double r, const std::vector<double>& strikes,
                                     const std::vector<double>& expiries,
                                     const std::vector<double>& call_prices,
                                     HestonParams init = {2.0, 0.04, 0.4, -0.5, 0.04}) {
    const std::size_t n = strikes.size();

    auto residual = [&](const std::vector<double>& q, std::vector<double>& out) {
        // Keep parameters in their valid domain via simple transforms.
        HestonParams hp{std::fabs(q[0]), std::fabs(q[1]), std::fabs(q[2]),
                        std::max(-0.999, std::min(0.999, q[3])), std::fabs(q[4])};
        for (std::size_t k = 0; k < n; ++k)
            out[k] = heston_call(hp, S, strikes[k], r, expiries[k]) - call_prices[k];
    };

    std::vector<double> p0 = {init.kappa, init.theta, init.sigma, init.rho, init.v0};
    const opt::LMResult res = opt::levenberg_marquardt(residual, p0, static_cast<int>(n));
    return {std::fabs(res.params[0]), std::fabs(res.params[1]), std::fabs(res.params[2]),
            std::max(-0.999, std::min(0.999, res.params[3])), std::fabs(res.params[4])};
}

}  // namespace pricer

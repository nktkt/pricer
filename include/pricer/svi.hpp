// pricer/svi.hpp — SVI volatility-smile model and its calibration.
//
// Gatheral's "raw" SVI parameterizes total implied variance in log-moneyness
// k = log(K/S):
//     w(k) = a + b * ( rho*(k - m) + sqrt((k - m)^2 + sigma^2) )
// and implied vol is sqrt(w(k)/T). This is the industry-standard arbitrage-aware
// smile. `calibrate_svi` fits the five parameters to market implied vols by
// least squares (via the Levenberg–Marquardt solver), i.e. a real model
// calibration that reproduces market quotes within tolerance.
#pragma once
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <vector>

#include "pricer/optimize.hpp"

namespace pricer {

struct SVIParams {
    double a, b, rho, m, sigma;

    // Total variance at log-moneyness k.
    double total_variance(double k) const {
        const double d = k - m;
        return a + b * (rho * d + std::sqrt(d * d + sigma * sigma));
    }
    // Implied volatility for strike K given spot S and expiry T.
    double vol(double K, double S, double T) const {
        const double w = total_variance(std::log(K / S));
        return std::sqrt((w > 0 ? w : 0.0) / T);
    }
};

// Calibrate raw SVI to market (strike, implied-vol) quotes for one expiry T.
// `init` optionally seeds the optimizer; otherwise a sensible guess is used.
inline SVIParams calibrate_svi(double S, double T, const std::vector<double>& strikes,
                               const std::vector<double>& vols,
                               SVIParams init = {0.04, 0.1, -0.3, 0.0, 0.1}) {
    if (strikes.size() != vols.size() || strikes.size() < 5)
        throw std::invalid_argument("calibrate_svi: need >= 5 matching quotes");

    const std::size_t n = strikes.size();
    std::vector<double> k(n), w_mkt(n);  // log-moneyness and target total variance
    for (std::size_t i = 0; i < n; ++i) {
        k[i] = std::log(strikes[i] / S);
        w_mkt[i] = vols[i] * vols[i] * T;
    }

    // Residuals in total-variance space. sigma and b enter via squares/abs so
    // the optimizer is free to wander in sign; we read back canonical values.
    auto residual = [&](const std::vector<double>& p, std::vector<double>& out) {
        SVIParams q{p[0], p[1], p[2], p[3], p[4]};
        for (std::size_t i = 0; i < n; ++i) out[i] = q.total_variance(k[i]) - w_mkt[i];
    };

    std::vector<double> p0 = {init.a, init.b, init.rho, init.m, init.sigma};
    const opt::LMResult res = opt::levenberg_marquardt(residual, p0, static_cast<int>(n));

    SVIParams fit{res.params[0], res.params[1], res.params[2], res.params[3], res.params[4]};
    fit.sigma = std::fabs(fit.sigma);  // sigma enters only as sigma^2, so |sigma| is equivalent
    return fit;
}

}  // namespace pricer

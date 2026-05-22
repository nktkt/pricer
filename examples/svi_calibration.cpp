// svi_calibration.cpp — calibrating a volatility smile to market quotes.
// We build a synthetic market from a known "true" SVI smile, add a tiny
// perturbation to each quote so the fit is a genuine least-squares problem,
// and then recover the smile with `calibrate_svi` (Levenberg-Marquardt).
#include "pricer/svi.hpp"
#include <cstdio>
#include <cmath>
#include <vector>

using namespace pricer;

int main() {
    // Market setup: spot and a single one-year expiry.
    const double S = 100.0, T = 1.0;

    // The "true" smile we pretend the market is generated from.
    const SVIParams truth{0.04, 0.10, -0.40, 0.00, 0.20};

    // Synthetic market: a strip of strikes around the money.
    const std::vector<double> strikes = {70, 80, 90, 100, 110, 120, 130, 140};

    // Build market implied vols: the true vol plus a small deterministic
    // perturbation (up to ~0.002 vol points) so calibration is a real fit,
    // not an exact interpolation.
    std::vector<double> market_vols(strikes.size());
    for (std::size_t i = 0; i < strikes.size(); ++i) {
        const double true_vol = truth.vol(strikes[i], S, T);
        const double bump = 0.002 * std::sin(static_cast<double>(i));  // fixed pattern
        market_vols[i] = true_vol + bump;
    }

    // Calibrate the five SVI parameters to the noisy quotes.
    const SVIParams fit = calibrate_svi(S, T, strikes, market_vols);

    // Compare fitted vols against the market quotes, strike by strike.
    std::printf("%8s | %12s | %12s | %12s\n", "K", "market vol", "fitted vol", "abs error");
    std::printf("---------|--------------|--------------|-------------\n");
    double sse = 0.0;
    for (std::size_t i = 0; i < strikes.size(); ++i) {
        const double fv = fit.vol(strikes[i], S, T);
        const double err = std::fabs(fv - market_vols[i]);
        sse += err * err;
        std::printf("%8.1f | %12.6f | %12.6f | %12.6f\n",
                    strikes[i], market_vols[i], fv, err);
    }
    const double rms = std::sqrt(sse / static_cast<double>(strikes.size()));

    // Report the recovered parameters and the overall fit quality.
    std::printf("\nfitted SVI parameters:\n");
    std::printf("  a     = %10.6f\n", fit.a);
    std::printf("  b     = %10.6f\n", fit.b);
    std::printf("  rho   = %10.6f\n", fit.rho);
    std::printf("  m     = %10.6f\n", fit.m);
    std::printf("  sigma = %10.6f\n", fit.sigma);
    std::printf("\nRMS fit error = %.6f vol points\n", rms);

    std::printf("\nSVI calibrated to the smile by least squares (Levenberg-Marquardt),\n");
    std::printf("reproducing the market quotes within a small RMS error.\n");
    return 0;
}

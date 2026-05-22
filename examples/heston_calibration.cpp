// heston_calibration.cpp — fitting Heston to a market option grid.
// Build a synthetic "market" of European call quotes from a known Heston model
// (plus a tiny perturbation), then recover the five parameters by least squares
// (Levenberg–Marquardt) and check how closely the fitted prices match the market.
#include "pricer/heston.hpp"
#include <cstdio>
#include <cmath>
#include <vector>

using namespace pricer;

int main() {
    const double S = 100.0, r = 0.05;

    // The "true" model that generated the market — calibration should recover it.
    const HestonParams truth{2.0, 0.04, 0.30, -0.60, 0.04};

    // Market grid: five strikes at two expiries → ten call quotes.
    const std::vector<double> strikes_grid = {80.0, 90.0, 100.0, 110.0, 120.0};
    const std::vector<double> expiries_grid = {0.5, 1.0};

    // Flatten the grid into three parallel vectors (length 10), in (T, K) order.
    std::vector<double> strikes_vec, expiries_vec, prices_vec;
    int index = 0;
    for (double T : expiries_grid) {
        for (double K : strikes_grid) {
            // Market price = true Heston price plus a small deterministic wobble,
            // so the fit is a genuine least-squares problem, not an exact solve.
            const double clean = heston_call(truth, S, K, r, T);
            const double market = clean * (1.0 + 0.0005 * std::sin(static_cast<double>(index)));
            strikes_vec.push_back(K);
            expiries_vec.push_back(T);
            prices_vec.push_back(market);
            ++index;
        }
    }

    std::printf("Calibrating Heston to a %zu-quote market grid (S=%.1f, r=%.2f)\n",
                prices_vec.size(), S, r);
    std::printf("True params:   kappa=%.4f theta=%.4f sigma=%.4f rho=%.4f v0=%.4f\n\n",
                truth.kappa, truth.theta, truth.sigma, truth.rho, truth.v0);

    // Run the calibration (Levenberg–Marquardt, default starting guess).
    const HestonParams fit =
        calibrate_heston(S, r, strikes_vec, expiries_vec, prices_vec);

    std::printf("Fitted params: kappa=%.4f theta=%.4f sigma=%.4f rho=%.4f v0=%.4f\n\n",
                fit.kappa, fit.theta, fit.sigma, fit.rho, fit.v0);

    // Compare market vs. fitted prices across the grid.
    std::printf("%8s | %6s | %12s | %12s | %12s\n", "K", "T", "market", "fitted", "abs err");
    std::printf("---------|--------|--------------|--------------|--------------\n");

    double sse = 0.0, max_abs = 0.0;
    for (std::size_t k = 0; k < prices_vec.size(); ++k) {
        const double model = heston_call(fit, S, strikes_vec[k], r, expiries_vec[k]);
        const double err = std::fabs(model - prices_vec[k]);
        sse += err * err;
        if (err > max_abs) max_abs = err;
        std::printf("%8.2f | %6.2f | %12.6f | %12.6f | %12.6f\n",
                    strikes_vec[k], expiries_vec[k], prices_vec[k], model, err);
    }
    const double rms = std::sqrt(sse / static_cast<double>(prices_vec.size()));

    std::printf("\nRMS price error = %.6f   max abs price error = %.6f\n", rms, max_abs);
    std::printf("\nHeston calibrated to the option grid by least squares "
                "(Levenberg-Marquardt).\n");
    return 0;
}

// sabr_smile.cpp — the SABR model: a whole volatility smile from one formula,
// calibration to market quotes, and a Monte Carlo cross-check.
//
// SABR (Hagan et al. 2002) gives the Black implied vol of an option on a forward
// in closed form, so a single set of parameters (α, β, ρ, ν) describes the entire
// smile. This demo prints the SABR smile, recovers the parameters from quotes by
// least-squares calibration, and confirms the asymptotic price against a direct
// simulation of the SABR SDE.
#include "pricer/sabr.hpp"

#include <cstdio>
#include <vector>

using namespace pricer;

int main() {
    const double F = 100.0, T = 1.0, df = 1.0;  // option on a forward
    // alpha, beta, rho, nu. With beta=0.5, alpha is a CEV-level vol; alpha=2 puts
    // the ATM lognormal vol near alpha/F^(1-beta) = 2/10 = 20%, a realistic smile.
    const SabrParams p{2.0, 0.5, -0.30, 0.40};

    std::printf("SABR smile   F=%.0f T=%.1f   alpha=%.2f beta=%.2f rho=%.2f nu=%.2f\n\n",
                F, T, p.alpha, p.beta, p.rho, p.nu);
    std::printf("  strike   implied vol   Black-76 call\n");
    const std::vector<double> strikes{70, 80, 90, 100, 110, 120, 140};
    for (double K : strikes) {
        const double iv = sabr_implied_vol(F, K, T, p);
        const double c = sabr_black_price(OptionType::Call, F, K, T, df, p);
        std::printf("  %6.0f   %9.4f     %9.4f\n", K, iv, c);
    }

    // --- Calibration: recover (alpha, rho, nu) from the smile above (fixed beta) ---
    std::vector<double> vols;
    for (double K : strikes) vols.push_back(sabr_implied_vol(F, K, T, p));
    const SabrFit fit = calibrate_sabr(F, T, strikes, vols, /*beta=*/0.5);
    std::printf("\nCalibration (beta fixed at 0.5):\n");
    std::printf("  recovered  alpha=%.4f rho=%.4f nu=%.4f   (RMS vol err %.2e)\n",
                fit.params.alpha, fit.params.rho, fit.params.nu, fit.rms_vol_error);

    // --- Monte Carlo of the SABR SDE vs the Hagan/Black-76 ATM price ---
    const double an = sabr_black_price(OptionType::Call, F, 100.0, T, df, p);
    const double mc = sabr_price_mc(OptionType::Call, F, 100.0, T, df, p, 200, 500'000);
    std::printf("\nATM call   Hagan/Black-76: %.4f   SABR-SDE Monte Carlo: %.4f   (diff %+.4f)\n",
                an, mc, mc - an);
    return 0;
}

// Tests for the SABR stochastic-volatility model (sabr.hpp): the Hagan implied-
// vol approximation, Black-76 pricing, least-squares calibration, and a Monte
// Carlo cross-check of the asymptotic formula against the simulated SDE.
#include "check.hpp"
#include "pricer/black_scholes.hpp"
#include "pricer/sabr.hpp"

using namespace pricer;

int main() {
    const double F = 100.0, T = 1.0, df = 1.0;  // work on the forward; df folded out

    // --- Reduction: beta = 1, nu = 0 gives a flat lognormal smile sigma = alpha ---
    const SabrParams flat{0.25, 1.0, 0.0, 0.0};
    for (double K : {80.0, 100.0, 125.0})
        check::approx("flat smile == alpha", sabr_implied_vol(F, K, T, flat), 0.25, 1e-10);

    // --- ATM branch is the continuous limit of the general branch as K -> F ---
    const SabrParams p{0.20, 0.5, -0.3, 0.4};
    const double atm = sabr_implied_vol(F, F, T, p);
    check::approx("general branch -> ATM as K->F", sabr_implied_vol(F, F * (1.0 + 1e-7), T, p), atm,
                  1e-5);

    // --- Skew: negative rho lifts low-strike vols above high-strike vols ---
    const SabrParams skew{0.20, 0.5, -0.5, 0.3};
    check::is_true("rho<0 gives downward skew",
                   sabr_implied_vol(F, 80.0, T, skew) > sabr_implied_vol(F, 120.0, T, skew));

    // --- Convexity: positive vol-of-vol lifts the wings above the ATM (a smile) ---
    const SabrParams smile{0.20, 0.5, 0.0, 0.6};
    const double atm_s = sabr_implied_vol(F, 100.0, T, smile);
    check::is_true("nu>0 gives a smile (wings up)",
                   sabr_implied_vol(F, 70.0, T, smile) > atm_s &&
                       sabr_implied_vol(F, 140.0, T, smile) > atm_s);

    // --- Calibration round-trip: recover (alpha, rho, nu) from model-generated vols ---
    const SabrParams truth{0.22, 0.5, -0.35, 0.45};
    const std::vector<double> strikes{70, 85, 100, 115, 130, 150};
    std::vector<double> vols;
    for (double K : strikes) vols.push_back(sabr_implied_vol(F, K, T, truth));
    const SabrFit fit = calibrate_sabr(F, T, strikes, vols, /*beta=*/0.5);
    check::approx("calibrated alpha", fit.params.alpha, truth.alpha, 1e-3);
    check::approx("calibrated rho", fit.params.rho, truth.rho, 1e-3);
    check::approx("calibrated nu", fit.params.nu, truth.nu, 1e-3);
    check::is_true("calibration RMS vol error tiny", fit.rms_vol_error < 1e-5);

    // --- Black-76 pricing: positive value and forward put-call parity ---
    const double call = sabr_black_price(OptionType::Call, F, 100.0, T, df, p);
    const double put = sabr_black_price(OptionType::Put, F, 100.0, T, df, p);
    check::is_true("ATM SABR call positive", call > 0.0);
    check::approx("forward put-call parity", call - put, df * (F - 100.0), 1e-9);

    // --- Monte Carlo of the SABR SDE cross-checks the Hagan/Black-76 price ---
    // (Hagan is an asymptotic approximation and the SDE is simulated, so this is a
    // moderate-tolerance sanity check, not an identity.)
    const SabrParams pm{0.20, 0.5, -0.2, 0.3};
    const double an_atm = sabr_black_price(OptionType::Call, F, 100.0, T, df, pm);
    const double mc_atm = sabr_price_mc(OptionType::Call, F, 100.0, T, df, pm, 200, 400'000);
    check::approx("SABR MC vs Hagan/Black76 (ATM)", mc_atm, an_atm, 0.25);

    // MC is reproducible for a fixed seed.
    const double mc_atm2 = sabr_price_mc(OptionType::Call, F, 100.0, T, df, pm, 200, 400'000);
    check::is_true("SABR MC reproducible", mc_atm == mc_atm2);

    return check::report("sabr");
}

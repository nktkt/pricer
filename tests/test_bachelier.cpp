// Tests for the Bachelier (normal) model (bachelier.hpp). The closed form is
// pinned by put-call parity, the at-the-money formula and the zero-vol intrinsic
// limit; the model is exercised on a negative-rate scenario (where Black-Scholes
// can't go); Greeks are checked against finite differences and the implied-vol
// solver round-trips.
#include "check.hpp"
#include "pricer/bachelier.hpp"
#include "pricer/black_scholes.hpp"

using namespace pricer;

int main() {
    const double F = 100, K = 100, sigma_n = 15.0, T = 1.0, df = 0.95;

    // --- Put-call parity: C - P = df * (F - K) ---
    const double c = bachelier_price(OptionType::Call, F, 105.0, sigma_n, T, df);
    const double p = bachelier_price(OptionType::Put, F, 105.0, sigma_n, T, df);
    check::approx("put-call parity", c - p, df * (F - 105.0), 1e-12);

    // --- At the money: price = df * sigma_n * sqrt(T) * phi(0) ---
    constexpr double phi0 = 0.3989422804014327;
    check::approx("ATM call closed form", bachelier_price(OptionType::Call, F, K, sigma_n, T, df),
                  df * sigma_n * std::sqrt(T) * phi0, 1e-12);
    check::approx("ATM call == ATM put", bachelier_price(OptionType::Call, F, K, sigma_n, T, df),
                  bachelier_price(OptionType::Put, F, K, sigma_n, T, df), 1e-12);

    // --- sigma_n -> 0 collapses to the discounted intrinsic on the forward ---
    check::approx("zero-vol call -> intrinsic", bachelier_price(OptionType::Call, 110.0, K, 1e-9, T, df),
                  df * 10.0, 1e-6);
    check::approx("zero-vol OTM call -> 0", bachelier_price(OptionType::Call, 90.0, K, 1e-9, T, df),
                  0.0, 1e-6);

    // --- Negative rates: F and K below zero, where the normal model is needed ---
    const double Fn = -0.005, Kn = 0.0, sn = 0.01;  // -0.5% forward rate, ATM-ish floor
    const double neg = bachelier_price(OptionType::Put, Fn, Kn, sn, T, 1.0);
    check::is_true("negative-rate put is positive", neg > 0.0);
    check::approx("negative-rate put parity",
                  bachelier_price(OptionType::Call, Fn, Kn, sn, T, 1.0) - neg, (Fn - Kn), 1e-12);

    // --- Implied vol round-trip recovers the input normal vol ---
    const double price = bachelier_price(OptionType::Call, F, 108.0, sigma_n, T, df);
    check::approx("implied normal vol round-trip",
                  bachelier_implied_vol(OptionType::Call, price, F, 108.0, T, df), sigma_n, 1e-8);
    // ATM round-trip too (the seed is exact there).
    const double pa = bachelier_price(OptionType::Put, F, K, sigma_n, T, df);
    check::approx("ATM implied vol round-trip",
                  bachelier_implied_vol(OptionType::Put, pa, F, K, T, df), sigma_n, 1e-8);

    // --- Greeks vs central finite differences ---
    const Greeks g = bachelier_greeks(OptionType::Call, F, 102.0, sigma_n, T, df);
    const double hF = 0.01, hs = 0.001, hT = 1e-5;
    const double dlt = (bachelier_price(OptionType::Call, F + hF, 102.0, sigma_n, T, df) -
                        bachelier_price(OptionType::Call, F - hF, 102.0, sigma_n, T, df)) / (2 * hF);
    check::approx("delta vs FD", g.delta, dlt, 1e-6);
    const double gam = (bachelier_price(OptionType::Call, F + hF, 102.0, sigma_n, T, df) -
                        2 * bachelier_price(OptionType::Call, F, 102.0, sigma_n, T, df) +
                        bachelier_price(OptionType::Call, F - hF, 102.0, sigma_n, T, df)) / (hF * hF);
    check::approx("gamma vs FD", g.gamma, gam, 1e-5);
    const double veg = (bachelier_price(OptionType::Call, F, 102.0, sigma_n + hs, T, df) -
                        bachelier_price(OptionType::Call, F, 102.0, sigma_n - hs, T, df)) / (2 * hs);
    check::approx("vega vs FD", g.vega, veg, 1e-6);
    const double tht = -(bachelier_price(OptionType::Call, F, 102.0, sigma_n, T + hT, df) -
                         bachelier_price(OptionType::Call, F, 102.0, sigma_n, T - hT, df)) / (2 * hT);
    check::approx("theta vs FD", g.theta, tht, 1e-4);

    // --- Monte Carlo matches the closed form, including the negative-rate case ---
    check::approx("MC call vs closed form",
                  bachelier_price_mc(OptionType::Call, F, 102.0, sigma_n, T, 2'000'000, 12345, df),
                  bachelier_price(OptionType::Call, F, 102.0, sigma_n, T, df), 0.05);
    check::approx("MC negative-rate put",
                  bachelier_price_mc(OptionType::Put, Fn, Kn, sn, T, 2'000'000, 7, 1.0), neg, 5e-5);

    // --- Reproducible for a fixed seed ---
    check::is_true("Bachelier MC reproducible",
                   bachelier_price_mc(OptionType::Call, F, K, sigma_n, T, 100000) ==
                       bachelier_price_mc(OptionType::Call, F, K, sigma_n, T, 100000));

    return check::report("bachelier");
}

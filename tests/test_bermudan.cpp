// Tests for Bermudan option pricing (bermudan.hpp). A Bermudan sits between a
// European and an American option, so the LSM price is checked against both
// limits: one exercise date reproduces the European price, many equally-spaced
// dates approach the American (binomial) price, and adding dates only raises the
// value.
#include "check.hpp"
#include "pricer/american.hpp"
#include "pricer/bermudan.hpp"
#include "pricer/black_scholes.hpp"

using namespace pricer;

int main() {
    const double S = 100, K = 100, r = 0.05, sigma = 0.20, T = 1.0;
    const long n = 200'000;

    const double eur_put = black_scholes_put(S, K, r, sigma, T);
    const double am_put = binomial_price(OptionType::Put, S, K, r, sigma, T, 2000, true);

    // --- A single exercise date (expiry only) reproduces the European price ---
    const double berm1 = bermudan_lsm(OptionType::Put, S, K, r, sigma, equally_spaced_dates(T, 1), n);
    check::approx("Bermudan(1 date) == European put", berm1, eur_put, 0.06);

    // --- A few exercise dates: European < Bermudan < American ---
    const double berm4 = bermudan_lsm(OptionType::Put, S, K, r, sigma, equally_spaced_dates(T, 4), n);
    check::is_true("Bermudan(4) > European put", berm4 > eur_put);
    check::is_true("Bermudan(4) < American put", berm4 < am_put + 0.05);

    // --- Adding exercise dates only raises the value (monotone toward American) ---
    const double berm12 = bermudan_lsm(OptionType::Put, S, K, r, sigma, equally_spaced_dates(T, 12), n);
    const double berm50 = bermudan_lsm(OptionType::Put, S, K, r, sigma, equally_spaced_dates(T, 50), n);
    check::is_true("Bermudan rises with more dates: 4<=12", berm12 > berm4 - 0.05);
    check::is_true("Bermudan rises with more dates: 12<=50", berm50 > berm12 - 0.05);

    // --- Many equally-spaced dates approach the American (binomial) price ---
    check::approx("Bermudan(50 dates) ~ American put", berm50, am_put, 0.15);
    check::is_true("Bermudan(50) <= American put (+noise)", berm50 < am_put + 0.1);

    // --- A Bermudan call on a non-dividend stock equals the European call ---
    // (early exercise is never optimal, so the extra dates add nothing).
    const double eur_call = black_scholes_call(S, K, r, sigma, T);
    const double berm_call = bermudan_lsm(OptionType::Call, S, K, r, sigma, equally_spaced_dates(T, 12), n);
    check::approx("Bermudan call == European call (no div)", berm_call, eur_call, 0.12);

    // --- An irregular (non-uniform) schedule works and lands between the limits ---
    const std::vector<double> irregular{0.1, 0.35, 0.6, 0.85, 1.0};
    const double berm_irr = bermudan_lsm(OptionType::Put, S, K, r, sigma, irregular, n);
    check::is_true("irregular-schedule Bermudan between European and American",
                   berm_irr > eur_put && berm_irr < am_put + 0.05);

    // --- Reproducible for a fixed seed ---
    const double berm4b = bermudan_lsm(OptionType::Put, S, K, r, sigma, equally_spaced_dates(T, 4), n);
    check::is_true("Bermudan reproducible run-to-run", berm4 == berm4b);

    return check::report("bermudan");
}

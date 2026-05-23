// Tests for American (early-exercise) option pricing (american.hpp).
#include "check.hpp"
#include "pricer/american.hpp"
#include "pricer/black_scholes.hpp"

using namespace pricer;

int main() {
    const double S = 100, K = 100, r = 0.05, sigma = 0.20, T = 1.0;

    // --- binomial European converges to Black–Scholes ---
    const double bin_eur_c = binomial_price(OptionType::Call, S, K, r, sigma, T, 1000, false);
    const double bin_eur_p = binomial_price(OptionType::Put, S, K, r, sigma, T, 1000, false);
    check::approx("binomial European call vs BS", bin_eur_c, black_scholes_call(S, K, r, sigma, T),
                  0.01);
    check::approx("binomial European put vs BS", bin_eur_p, black_scholes_put(S, K, r, sigma, T),
                  0.01);

    // --- American call on a non-dividend stock equals the European call ---
    // (early exercise is never optimal, so the prices coincide).
    const double bin_am_c = binomial_price(OptionType::Call, S, K, r, sigma, T, 1000, true);
    check::approx("American call == European call (no div)", bin_am_c, bin_eur_c, 1e-9);

    // --- American put carries a positive early-exercise premium ---
    const double bin_am_p = binomial_price(OptionType::Put, S, K, r, sigma, T, 1000, true);
    check::is_true("American put > European put", bin_am_p > bin_eur_p);
    check::is_true("early-exercise premium is material", bin_am_p - bin_eur_p > 0.1);

    // --- a deep in-the-money American put is worth ~ its intrinsic value ---
    // (immediate exercise is optimal), and never less than intrinsic.
    const double deep = binomial_price(OptionType::Put, 60.0, K, r, sigma, T, 1000, true);
    check::is_true("deep-ITM American put >= intrinsic", deep >= 40.0 - 1e-6);
    check::approx("deep-ITM American put ~ intrinsic", deep, 40.0, 1.0);

    // --- Longstaff–Schwartz LSM agrees with the binomial American put ---
    const double lsm_p = lsm_american(OptionType::Put, S, K, r, sigma, T, 80'000, 50);
    check::approx("LSM American put vs binomial", lsm_p, bin_am_p, 0.1);

    // --- LSM American call (no div) matches the European call ---
    const double lsm_c = lsm_american(OptionType::Call, S, K, r, sigma, T, 80'000, 50);
    check::approx("LSM American call vs European", lsm_c, black_scholes_call(S, K, r, sigma, T), 0.15);

    // --- LSM is reproducible for a fixed seed ---
    const double lsm_p2 = lsm_american(OptionType::Put, S, K, r, sigma, T, 80'000, 50);
    check::is_true("LSM reproducible run-to-run", lsm_p == lsm_p2);

    // --- binomial Greeks (European) match the closed-form Greeks ---
    const Greeks bg = binomial_greeks(OptionType::Call, S, K, r, sigma, T, 2000, false);
    const Greeks cf = black_scholes_greeks(OptionType::Call, S, K, r, sigma, T);
    check::approx("binomial delta vs BS", bg.delta, cf.delta, 5e-3);
    check::approx("binomial gamma vs BS", bg.gamma, cf.gamma, 5e-3);
    check::approx("binomial vega vs BS", bg.vega, cf.vega, 0.2);
    check::approx("binomial theta vs BS", bg.theta, cf.theta, 0.1);
    check::approx("binomial rho vs BS", bg.rho, cf.rho, 0.5);

    // --- American Greeks: call (no div) matches European; put has the right signs ---
    const Greeks am_c = binomial_greeks(OptionType::Call, S, K, r, sigma, T, 2000, true);
    check::approx("American call delta == European", am_c.delta, bg.delta, 5e-3);
    const Greeks am_p = binomial_greeks(OptionType::Put, S, K, r, sigma, T, 2000, true);
    check::is_true("American put delta < 0", am_p.delta < 0.0);
    check::is_true("American put gamma > 0", am_p.gamma > 0.0);
    check::is_true("American put vega > 0", am_p.vega > 0.0);

    // --- binomial converges: more steps shrink the error vs Black–Scholes ---
    const double coarse = binomial_price(OptionType::Call, S, K, r, sigma, T, 50, false);
    const double fine = binomial_price(OptionType::Call, S, K, r, sigma, T, 2000, false);
    const double exact = black_scholes_call(S, K, r, sigma, T);
    check::is_true("binomial converges with steps",
                   std::fabs(fine - exact) < std::fabs(coarse - exact));

    return check::report("american");
}

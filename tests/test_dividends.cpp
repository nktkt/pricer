// Tests for continuous dividend yield (q) across the pricing stack.
#include "check.hpp"
#include "pricer/adjoint.hpp"
#include "pricer/american.hpp"
#include "pricer/black_scholes.hpp"
#include "pricer/greeks_ad.hpp"
#include "pricer/implied_vol.hpp"

#include <cmath>

using namespace pricer;

int main() {
    const double S = 100, K = 100, r = 0.05, sigma = 0.20, T = 1.0, q = 0.03;

    // --- dividend-adjusted put–call parity: C - P = S e^{-qT} - K e^{-rT} ---
    const double c = black_scholes_call(S, K, r, sigma, T, q);
    const double p = black_scholes_put(S, K, r, sigma, T, q);
    check::approx("dividend put-call parity", c - p,
                  S * std::exp(-q * T) - K * std::exp(-r * T), 1e-9);

    // --- a dividend lowers the call and raises the put; q=0 is plain BS ---
    check::is_true("dividend lowers the call", c < black_scholes_call(S, K, r, sigma, T, 0.0));
    check::is_true("dividend raises the put", p > black_scholes_put(S, K, r, sigma, T, 0.0));
    check::approx("q=0 reduces to plain BS", black_scholes_call(S, K, r, sigma, T, 0.0),
                  black_scholes_call(S, K, r, sigma, T), 1e-12);

    // --- AD and AAD Greeks with q match the closed-form Greeks with q ---
    const Greeks cf = black_scholes_greeks(OptionType::Call, S, K, r, sigma, T, q);
    const Greeks ad = black_scholes_greeks_ad(OptionType::Call, S, K, r, sigma, T, q);
    const Greeks aad = black_scholes_greeks_aad(OptionType::Call, S, K, r, sigma, T, q);
    check::approx("AD delta with q", ad.delta, cf.delta, 1e-7);
    check::approx("AD vega with q", ad.vega, cf.vega, 1e-6);
    check::approx("AD gamma with q", ad.gamma, cf.gamma, 1e-7);
    check::approx("AD theta with q", ad.theta, cf.theta, 1e-6);
    check::approx("AAD price with q", aad.price, cf.price, 1e-9);
    check::approx("AAD delta with q", aad.delta, cf.delta, 1e-7);
    check::approx("AAD rho with q", aad.rho, cf.rho, 1e-6);

    // --- implied vol round-trips through the dividend-aware price ---
    const double iv = implied_vol(OptionType::Call, c, S, K, r, T, 1e-10, 100, q);
    check::approx("implied vol round-trip with q", iv, sigma, 1e-6);

    // --- binomial European with q converges to the closed-form price ---
    const double bin_eur = binomial_price(OptionType::Call, S, K, r, sigma, T, 2000, false, q);
    check::approx("binomial European with q vs BS", bin_eur, c, 0.01);

    // --- dividends make American CALL early exercise valuable (premium > 0) ---
    const double qc = 0.10;  // a high yield makes early exercise of the call optimal
    const double am_call = binomial_price(OptionType::Call, S, K, r, sigma, T, 1000, true, qc);
    const double eu_call = binomial_price(OptionType::Call, S, K, r, sigma, T, 1000, false, qc);
    check::is_true("dividend American call premium > 0", am_call - eu_call > 0.01);

    // --- LSM American put with q agrees with the binomial American put with q ---
    const double bin_am_p = binomial_price(OptionType::Put, S, K, r, sigma, T, 1000, true, q);
    const double lsm_am_p = lsm_american(OptionType::Put, S, K, r, sigma, T, 80'000, 50, 12345, q);
    check::approx("LSM American put with q vs binomial", lsm_am_p, bin_am_p, 0.1);

    return check::report("dividends");
}

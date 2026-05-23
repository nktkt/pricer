// Tests for digital (binary) options (digital.hpp). The closed forms are pinned
// by the exact relationships that tie digitals to vanillas: the vanilla
// decomposition, cash/asset parities, and the digital-as-strike-derivative
// identity. The Monte Carlo engine cross-checks the closed form.
#include "check.hpp"
#include "pricer/black_scholes.hpp"
#include "pricer/digital.hpp"

using namespace pricer;

int main() {
    const double S = 100, K = 100, r = 0.05, sigma = 0.20, T = 1.0;

    // Shorthands.
    auto cash = [&](OptionType t, double q = 0.0) {
        return digital_price(t, DigitalType::CashOrNothing, S, K, r, sigma, T, 1.0, q);
    };
    auto asset = [&](OptionType t, double q = 0.0) {
        return digital_price(t, DigitalType::AssetOrNothing, S, K, r, sigma, T, 1.0, q);
    };

    // --- Vanilla decomposition: call = asset-or-nothing − K · cash-or-nothing(1) ---
    check::approx("call == asset_on - K*cash_on", asset(OptionType::Call) - K * cash(OptionType::Call),
                  black_scholes_call(S, K, r, sigma, T), 1e-10);
    check::approx("put == K*cash_on - asset_on", K * cash(OptionType::Put) - asset(OptionType::Put),
                  black_scholes_put(S, K, r, sigma, T), 1e-10);

    // --- Cash-or-nothing parity: a $1 bet pays out either way, worth e^{-rT} ---
    check::approx("cash call + cash put == e^{-rT}", cash(OptionType::Call) + cash(OptionType::Put),
                  std::exp(-r * T), 1e-12);
    // --- Asset-or-nothing parity: you always receive the asset, worth S·e^{-qT} ---
    check::approx("asset call + asset put == S", asset(OptionType::Call) + asset(OptionType::Put), S,
                  1e-10);
    const double qd = 0.03;
    check::approx("asset parity with q == S·e^{-qT}",
                  asset(OptionType::Call, qd) + asset(OptionType::Put, qd), S * std::exp(-qd * T),
                  1e-10);

    // --- A cash-or-nothing call (cash=1) is −∂C/∂K of the vanilla call ---
    const double h = 1e-3;
    const double dCdK = (black_scholes_call(S, K - h, r, sigma, T) -
                         black_scholes_call(S, K + h, r, sigma, T)) / (2 * h);
    check::approx("cash-or-nothing call == -dC/dK", cash(OptionType::Call), dCdK, 1e-5);

    // --- Vanilla decomposition still holds with a dividend yield q ---
    check::approx("call == asset_on - K*cash_on (q)",
                  asset(OptionType::Call, qd) - K * cash(OptionType::Call, qd),
                  black_scholes_call(S, K, r, sigma, T, qd), 1e-10);

    // --- Monte Carlo matches the closed form for every kind / type ---
    const long n = 4'000'000;
    check::approx("cash call MC", digital_price_mc(OptionType::Call, DigitalType::CashOrNothing, S, K,
                                                   r, sigma, T, n),
                  cash(OptionType::Call), 0.01);
    check::approx("cash put MC", digital_price_mc(OptionType::Put, DigitalType::CashOrNothing, S, K, r,
                                                  sigma, T, n),
                  cash(OptionType::Put), 0.01);
    check::approx("asset call MC", digital_price_mc(OptionType::Call, DigitalType::AssetOrNothing, S,
                                                    K, r, sigma, T, n),
                  asset(OptionType::Call), 0.2);
    check::approx("asset put MC", digital_price_mc(OptionType::Put, DigitalType::AssetOrNothing, S, K,
                                                   r, sigma, T, n),
                  asset(OptionType::Put), 0.2);

    // --- A cash-or-nothing call pays `cash` proportionally ---
    check::approx("cash scales linearly",
                  digital_price(OptionType::Call, DigitalType::CashOrNothing, S, K, r, sigma, T, 10.0),
                  10.0 * cash(OptionType::Call), 1e-12);

    // --- Reproducible for a fixed seed ---
    check::is_true("digital MC reproducible",
                   digital_price_mc(OptionType::Call, DigitalType::CashOrNothing, S, K, r, sigma, T,
                                    100000) ==
                       digital_price_mc(OptionType::Call, DigitalType::CashOrNothing, S, K, r, sigma,
                                        T, 100000));

    return check::report("digital");
}

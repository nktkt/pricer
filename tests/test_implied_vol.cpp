// Tests for the implied-volatility solver (round-trip and robustness).
#include "check.hpp"
#include "pricer/black_scholes.hpp"
#include "pricer/implied_vol.hpp"
#include <cstdio>
#include <initializer_list>

using namespace pricer;

int main() {
    const double S = 100, r = 0.05, T = 1.0;

    // Round-trip: price at a known vol, then recover that vol.
    for (double K : {80.0, 100.0, 120.0}) {
        for (double sig : {0.10, 0.20, 0.40}) {
            const double price = black_scholes_call(S, K, r, sig, T);
            const double iv = implied_vol(OptionType::Call, price, S, K, r, T);
            char nm[48];
            std::snprintf(nm, sizeof nm, "iv K=%.0f sig=%.2f", K, sig);
            check::approx(nm, iv, sig, 1e-6);
        }
    }

    // Puts too.
    const double pput = black_scholes_put(S, 100, r, 0.25, T);
    check::approx("iv put", implied_vol(OptionType::Put, pput, S, 100, r, T), 0.25, 1e-6);

    return check::report("implied_vol");
}

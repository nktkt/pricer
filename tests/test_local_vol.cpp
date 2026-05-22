// Tests for Dupire local volatility.
#include "check.hpp"
#include "pricer/black_scholes.hpp"
#include "pricer/local_vol.hpp"
#include <cstdio>
#include <initializer_list>

using namespace pricer;

int main() {
    const double S = 100, r = 0.05;

    // A FLAT implied-vol surface: prices generated at constant 0.20 must yield a
    // local vol of 0.20 everywhere (the classic Dupire sanity check).
    auto flat = [&](double K, double T) { return black_scholes_call(S, K, r, 0.20, T); };
    for (double K : {90.0, 100.0, 110.0}) {
        for (double T : {0.5, 1.0, 2.0}) {
            char nm[40];
            std::snprintf(nm, sizeof nm, "flat local vol K=%.0f T=%.1f", K, T);
            check::approx(nm, dupire_local_vol(flat, K, T, r), 0.20, 2e-3);
        }
    }

    return check::report("local_vol");
}

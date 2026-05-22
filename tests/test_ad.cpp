// Tests for forward-mode AD Greeks: they must equal the closed-form Greeks.
#include "check.hpp"
#include "pricer/black_scholes.hpp"
#include "pricer/dual.hpp"
#include "pricer/greeks_ad.hpp"
#include <cmath>
#include <cstdio>
#include <initializer_list>

using namespace pricer;

int main() {
    // Sanity on the dual core: d/dx (x^2 * exp(x)) at x=1.5 = (2x + x^2) e^x.
    {
        Dual<1, double> x = ad_var<1>(1.5, 0);
        Dual<1, double> f = x * x * exp(x);
        const double expect = (2 * 1.5 + 1.5 * 1.5) * std::exp(1.5);
        check::approx("dual derivative", f.d[0], expect, 1e-10);
    }

    // AD Greeks vs. closed form for several configurations.
    for (OptionType t : {OptionType::Call, OptionType::Put}) {
        const char* tn = (t == OptionType::Call) ? "call" : "put";
        for (double sig : {0.15, 0.30}) {
            const double S = 100, K = 105, r = 0.04, T = 0.75;
            const Greeks ad = black_scholes_greeks_ad(t, S, K, r, sig, T);
            const Greeks cf = black_scholes_greeks(t, S, K, r, sig, T);
            char nm[48];
            std::snprintf(nm, sizeof nm, "%s sig=%.2f price", tn, sig); check::approx(nm, ad.price, cf.price, 1e-10);
            std::snprintf(nm, sizeof nm, "%s sig=%.2f delta", tn, sig); check::approx(nm, ad.delta, cf.delta, 1e-10);
            std::snprintf(nm, sizeof nm, "%s sig=%.2f gamma", tn, sig); check::approx(nm, ad.gamma, cf.gamma, 1e-10);
            std::snprintf(nm, sizeof nm, "%s sig=%.2f vega",  tn, sig); check::approx(nm, ad.vega,  cf.vega,  1e-8);
            std::snprintf(nm, sizeof nm, "%s sig=%.2f theta", tn, sig); check::approx(nm, ad.theta, cf.theta, 1e-8);
            std::snprintf(nm, sizeof nm, "%s sig=%.2f rho",   tn, sig); check::approx(nm, ad.rho,   cf.rho,   1e-8);
        }
    }

    return check::report("ad");
}

// Tests for reverse-mode AAD Greeks against the closed form.
#include "check.hpp"
#include "pricer/adjoint.hpp"
#include "pricer/black_scholes.hpp"
#include <cstdio>
#include <initializer_list>

using namespace pricer;

int main() {
    // Tape sanity: f = (x*y + exp(x)) / y has known partials.
    {
        Tape t;
        Var x = make_var(t, 1.5), y = make_var(t, 2.0);
        Var f = (x * y + exp(x)) / y;
        auto adj = t.grad(f.idx);
        const double ex = std::exp(1.5);
        // df/dx = (y + e^x)/y ; df/dy = -(e^x)/y^2  (the x*y/y term has 0 net dy)
        check::approx("df/dx", adj[x.idx], (2.0 + ex) / 2.0, 1e-10);
        check::approx("df/dy", adj[y.idx], -ex / 4.0, 1e-10);
    }

    // AAD Greeks vs. closed form.
    for (OptionType t : {OptionType::Call, OptionType::Put}) {
        const char* tn = (t == OptionType::Call) ? "call" : "put";
        for (double sig : {0.15, 0.35}) {
            const double S = 100, K = 95, r = 0.03, T = 1.5;
            const Greeks ad = black_scholes_greeks_aad(t, S, K, r, sig, T);
            const Greeks cf = black_scholes_greeks(t, S, K, r, sig, T);
            char nm[48];
            std::snprintf(nm, sizeof nm, "%s sig=%.2f price", tn, sig); check::approx(nm, ad.price, cf.price, 1e-10);
            std::snprintf(nm, sizeof nm, "%s sig=%.2f delta", tn, sig); check::approx(nm, ad.delta, cf.delta, 1e-10);
            std::snprintf(nm, sizeof nm, "%s sig=%.2f vega",  tn, sig); check::approx(nm, ad.vega,  cf.vega,  1e-8);
            std::snprintf(nm, sizeof nm, "%s sig=%.2f theta", tn, sig); check::approx(nm, ad.theta, cf.theta, 1e-8);
            std::snprintf(nm, sizeof nm, "%s sig=%.2f rho",   tn, sig); check::approx(nm, ad.rho,   cf.rho,   1e-8);
            std::snprintf(nm, sizeof nm, "%s sig=%.2f gamma", tn, sig); check::approx(nm, ad.gamma, cf.gamma, 1e-10);
        }
    }

    return check::report("adjoint");
}

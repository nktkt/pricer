// Tests for closed-form Greeks: cross-check against central finite differences.
#include "check.hpp"
#include "pricer/black_scholes.hpp"
#include <cstdio>
#include <initializer_list>

using namespace pricer;

static double price(OptionType t, double S, double K, double r, double sig, double T) {
    return black_scholes_price(t, S, K, r, sig, T);
}

int main() {
    const double S = 100, K = 100, r = 0.05, sigma = 0.20, T = 1.0;

    for (OptionType t : {OptionType::Call, OptionType::Put}) {
        const char* nm = (t == OptionType::Call) ? "call" : "put";
        const Greeks g = black_scholes_greeks(t, S, K, r, sigma, T);

        const double hS = 0.01, hSig = 1e-4, hT = 1e-4, hr = 1e-4;
        const double fd_delta = (price(t, S + hS, K, r, sigma, T) - price(t, S - hS, K, r, sigma, T)) / (2 * hS);
        const double fd_gamma = (price(t, S + hS, K, r, sigma, T) - 2 * g.price + price(t, S - hS, K, r, sigma, T)) / (hS * hS);
        const double fd_vega  = (price(t, S, K, r, sigma + hSig, T) - price(t, S, K, r, sigma - hSig, T)) / (2 * hSig);
        const double fd_theta = -(price(t, S, K, r, sigma, T + hT) - price(t, S, K, r, sigma, T - hT)) / (2 * hT);
        const double fd_rho   = (price(t, S, K, r + hr, sigma, T) - price(t, S, K, r - hr, sigma, T)) / (2 * hr);

        char label[64];
        std::snprintf(label, sizeof label, "%s delta", nm); check::approx(label, g.delta, fd_delta, 1e-4);
        std::snprintf(label, sizeof label, "%s gamma", nm); check::approx(label, g.gamma, fd_gamma, 1e-4);
        std::snprintf(label, sizeof label, "%s vega",  nm); check::approx(label, g.vega,  fd_vega,  1e-2);
        std::snprintf(label, sizeof label, "%s theta", nm); check::approx(label, g.theta, fd_theta, 1e-2);
        std::snprintf(label, sizeof label, "%s rho",   nm); check::approx(label, g.rho,   fd_rho,   1e-2);
    }

    return check::report("greeks");
}

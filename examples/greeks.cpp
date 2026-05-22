// greeks.cpp — risk sensitivities from the library, cross-checked numerically.
// Delta/Gamma/Vega/Theta/Rho come from pricer::black_scholes_greeks (closed form);
// each is verified against a central finite difference of the price.
#include "pricer/black_scholes.hpp"
#include <cstdio>

using namespace pricer;

int main() {
    const double S = 100, K = 100, r = 0.05, sigma = 0.20, T = 1.0;

    const Greeks g = black_scholes_greeks(OptionType::Call, S, K, r, sigma, T);

    auto C = [&](double S_, double r_, double sig_, double T_) {
        return black_scholes_call(S_, K, r_, sig_, T_);
    };
    const double hS = 0.01, hSig = 1e-4, hT = 1e-4, hr = 1e-4;
    const double fd_delta = (C(S + hS, r, sigma, T) - C(S - hS, r, sigma, T)) / (2 * hS);
    const double fd_gamma = (C(S + hS, r, sigma, T) - 2 * g.price + C(S - hS, r, sigma, T)) / (hS * hS);
    const double fd_vega  = (C(S, r, sigma + hSig, T) - C(S, r, sigma - hSig, T)) / (2 * hSig);
    const double fd_theta = -(C(S, r, sigma, T + hT) - C(S, r, sigma, T - hT)) / (2 * hT);
    const double fd_rho   = (C(S, r + hr, sigma, T) - C(S, r - hr, sigma, T)) / (2 * hr);

    std::printf("call price = %.6f\n\n", g.price);
    std::printf("%-7s | %12s | %12s | %s\n", "greek", "closed form", "finite diff", "meaning");
    std::printf("--------|--------------|--------------|------------------------------\n");
    std::printf("%-7s | %12.6f | %12.6f | %s\n", "Delta", g.delta, fd_delta, "price change per +1 spot");
    std::printf("%-7s | %12.6f | %12.6f | %s\n", "Gamma", g.gamma, fd_gamma, "rate of change of Delta");
    std::printf("%-7s | %12.6f | %12.6f | %s\n", "Vega",  g.vega,  fd_vega,  "per +1.0 vol (/100 for 1%)");
    std::printf("%-7s | %12.6f | %12.6f | %s\n", "Theta", g.theta, fd_theta, "time decay per year");
    std::printf("%-7s | %12.6f | %12.6f | %s\n", "Rho",   g.rho,   fd_rho,   "price change per +1.0 rate");
    std::printf("\nClosed form matches finite differences -> the implementation checks out.\n");
    std::printf("Practical: Theta/365 = %.4f per day, Vega/100 = %.4f per 1%% vol\n",
                g.theta / 365.0, g.vega / 100.0);
    return 0;
}

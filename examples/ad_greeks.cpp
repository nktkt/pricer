// ad_greeks.cpp — exact Greeks by forward-mode automatic differentiation.
//
// One evaluation of Black–Scholes in dual numbers yields delta, vega, rho and
// theta together; a nested dual gives gamma. The results match the closed-form
// Greeks to machine precision — no bumping, no finite-difference error (contrast
// with examples/mc_greeks.cpp, which bumps a Monte Carlo price).
#include "pricer/black_scholes.hpp"
#include "pricer/greeks_ad.hpp"
#include <cmath>
#include <cstdio>

using namespace pricer;

int main() {
    const double S = 100, K = 100, r = 0.05, sigma = 0.20, T = 1.0;

    const Greeks ad = black_scholes_greeks_ad(OptionType::Call, S, K, r, sigma, T);
    const Greeks cf = black_scholes_greeks(OptionType::Call, S, K, r, sigma, T);

    std::printf("call: S=%.0f K=%.0f r=%.2f sigma=%.2f T=%.0f\n\n", S, K, r, sigma, T);
    std::printf("%-7s | %16s | %16s | %s\n", "greek", "auto-diff", "closed form", "abs diff");
    std::printf("--------|------------------|------------------|---------\n");
    auto row = [](const char* nm, double a, double c) {
        std::printf("%-7s | %16.10f | %16.10f | %.2e\n", nm, a, c, std::fabs(a - c));
    };
    row("price", ad.price, cf.price);
    row("delta", ad.delta, cf.delta);
    row("gamma", ad.gamma, cf.gamma);
    row("vega",  ad.vega,  cf.vega);
    row("theta", ad.theta, cf.theta);
    row("rho",   ad.rho,   cf.rho);
    std::printf("\nForward-mode AD reproduces every Greek to machine precision.\n");
    return 0;
}

// aad_greeks.cpp — all first-order Greeks in one reverse-mode (adjoint) pass.
//
// Reverse-mode AD records the price computation on a tape, then a single backward
// sweep yields delta, vega, rho and theta together — independent of the number
// of inputs. This is the production-grade complement to the forward-mode AD in
// examples/ad_greeks.cpp (gamma here still comes from forward AD).
#include "pricer/adjoint.hpp"
#include "pricer/black_scholes.hpp"
#include <cmath>
#include <cstdio>

using namespace pricer;

int main() {
    const double S = 100, K = 100, r = 0.05, sigma = 0.20, T = 1.0;

    const Greeks ad = black_scholes_greeks_aad(OptionType::Call, S, K, r, sigma, T);
    const Greeks cf = black_scholes_greeks(OptionType::Call, S, K, r, sigma, T);

    std::printf("call: S=%.0f K=%.0f r=%.2f sigma=%.2f T=%.0f\n", S, K, r, sigma, T);
    std::printf("(delta/vega/rho/theta from ONE adjoint backward sweep; gamma via forward AD)\n\n");
    std::printf("%-7s | %16s | %16s | %s\n", "greek", "adjoint AD", "closed form", "abs diff");
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
    std::printf("\nOne forward evaluation + one backward sweep gives every first-order Greek.\n");
    return 0;
}

// rainbow_options.cpp — two-asset rainbow (best-of / worst-of) options.
//
// A rainbow option pays on the better or worse performer of two correlated
// assets. Stulz (1982) prices these in closed form via the bivariate normal CDF.
// This demo prints best-of (max) and worst-of (min) call prices, checks the
// parity that ties them to two vanilla calls, and cross-checks against a
// correlated two-asset Monte Carlo. It also shows how correlation moves the prices.
#include "pricer/black_scholes.hpp"
#include "pricer/rainbow.hpp"

#include <cstdio>

using namespace pricer;

int main() {
    const double S1 = 100, S2 = 95, K = 100, r = 0.05, sig1 = 0.20, sig2 = 0.30, rho = 0.4, T = 1.0;
    const long n = 1'000'000;

    std::printf("Rainbow options on 2 assets   S1=%.0f S2=%.0f K=%.0f sigma1=%.2f sigma2=%.2f "
                "rho=%.1f T=%.1f\n\n",
                S1, S2, K, sig1, sig2, rho, T);

    const double cmax = rainbow_price(OptionType::Call, RainbowType::Max, S1, S2, K, r, sig1, sig2, rho, T);
    const double cmin = rainbow_price(OptionType::Call, RainbowType::Min, S1, S2, K, r, sig1, sig2, rho, T);
    const double cmax_mc = rainbow_price_mc(OptionType::Call, RainbowType::Max, S1, S2, K, r, sig1, sig2, rho, T, n);
    const double cmin_mc = rainbow_price_mc(OptionType::Call, RainbowType::Min, S1, S2, K, r, sig1, sig2, rho, T, n);
    const double c1 = black_scholes_call(S1, K, r, sig1, T), c2 = black_scholes_call(S2, K, r, sig2, T);

    std::printf("  best-of  (max) call : closed %.4f   MC %.4f\n", cmax, cmax_mc);
    std::printf("  worst-of (min) call : closed %.4f   MC %.4f\n", cmin, cmin_mc);
    std::printf("  vanilla calls       : c(S1)=%.4f  c(S2)=%.4f\n", c1, c2);
    std::printf("  parity  Cmax + Cmin = c1 + c2 :  %.4f = %.4f\n\n", cmax + cmin, c1 + c2);

    std::printf("  correlation sweep (best-of / worst-of call):\n");
    for (double rr : {-0.5, 0.0, 0.5, 0.9}) {
        const double a = rainbow_price(OptionType::Call, RainbowType::Max, S1, S2, K, r, sig1, sig2, rr, T);
        const double b = rainbow_price(OptionType::Call, RainbowType::Min, S1, S2, K, r, sig1, sig2, rr, T);
        std::printf("    rho=%+.1f : best-of %.4f   worst-of %.4f\n", rr, a, b);
    }
    return 0;
}

// bachelier_demo.cpp — the Bachelier (normal) model, and why it exists.
//
// Black–Scholes assumes the underlying is lognormal, so it can never price an
// option when the forward or strike is negative — but interest rates and spreads
// routinely are. The Bachelier model puts arithmetic Brownian motion on the
// forward (F_T = F + sigma_n*sqrt(T)*Z), which handles negative values cleanly.
// This demo prices a normal-vol option (closed form vs. Monte Carlo), recovers
// the normal vol by inversion, and prices a negative-rate floorlet.
#include "pricer/bachelier.hpp"
#include "pricer/black_scholes.hpp"

#include <cstdio>

using namespace pricer;

int main() {
    const double F = 100, K = 100, sigma_n = 15.0, T = 1.0, df = 0.95;

    std::printf("Bachelier (normal) model   F=%.0f K=%.0f sigma_n=%.1f T=%.1f df=%.2f\n\n",
                F, K, sigma_n, T, df);

    const double c = bachelier_price(OptionType::Call, F, K, sigma_n, T, df);
    const double c_mc = bachelier_price_mc(OptionType::Call, F, K, sigma_n, T, 4'000'000, 12345, df);
    const Greeks g = bachelier_greeks(OptionType::Call, F, K, sigma_n, T, df);
    std::printf("ATM call\n");
    std::printf("  closed form : %.4f   (= df*sigma_n*sqrt(T)*phi(0))\n", c);
    std::printf("  Monte Carlo : %.4f\n", c_mc);
    std::printf("  Greeks      : delta=%.4f gamma=%.5f vega=%.4f theta=%.4f\n",
                g.delta, g.gamma, g.vega, g.theta);
    std::printf("  implied normal vol from price: %.4f\n\n",
                bachelier_implied_vol(OptionType::Call, c, F, K, T, df));

    // --- A negative-rate floorlet: a put on a forward rate that is below zero ---
    const double Fr = -0.004, Kr = 0.0, sr = 0.012;  // -0.40% forward rate, 0% strike
    std::printf("Negative-rate floorlet (Black-Scholes cannot price this)\n");
    std::printf("  forward rate F=%.3f%%  strike K=%.1f%%  normal vol=%.2f%%\n",
                Fr * 100, Kr * 100, sr * 100);
    const double floor = bachelier_price(OptionType::Put, Fr, Kr, sr, T, 1.0);
    const double floor_mc = bachelier_price_mc(OptionType::Put, Fr, Kr, sr, T, 4'000'000, 7, 1.0);
    std::printf("  floorlet value: closed %.6f   MC %.6f\n", floor, floor_mc);
    return 0;
}

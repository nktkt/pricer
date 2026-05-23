// digital_options.cpp — digital (binary) options and how they build a vanilla.
//
// A digital option pays a discontinuous amount depending only on whether it
// finishes in the money: cash-or-nothing pays a fixed sum, asset-or-nothing pays
// the underlying. This demo prices both (closed form vs. Monte Carlo) and shows
// the exact decomposition vanilla call = asset-or-nothing − K · cash-or-nothing.
#include "pricer/black_scholes.hpp"
#include "pricer/digital.hpp"

#include <cstdio>

using namespace pricer;

int main() {
    const double S = 100, K = 100, r = 0.05, sigma = 0.20, T = 1.0;
    const long n = 4'000'000;

    std::printf("Digital options   S=%.0f K=%.0f r=%.2f sigma=%.2f T=%.1f   (%ld MC paths)\n\n",
                S, K, r, sigma, T, n);

    const double cash_c = digital_price(OptionType::Call, DigitalType::CashOrNothing, S, K, r, sigma, T);
    const double cash_c_mc = digital_price_mc(OptionType::Call, DigitalType::CashOrNothing, S, K, r, sigma, T, n);
    const double asset_c = digital_price(OptionType::Call, DigitalType::AssetOrNothing, S, K, r, sigma, T);
    const double asset_c_mc = digital_price_mc(OptionType::Call, DigitalType::AssetOrNothing, S, K, r, sigma, T, n);

    std::printf("CASH-OR-NOTHING call ($1 if S_T > K)\n");
    std::printf("  closed form : %.4f   (= e^-rT * P[S_T > K])\n", cash_c);
    std::printf("  Monte Carlo : %.4f\n\n", cash_c_mc);
    std::printf("ASSET-OR-NOTHING call (S_T if S_T > K)\n");
    std::printf("  closed form : %.4f\n", asset_c);
    std::printf("  Monte Carlo : %.4f\n\n", asset_c_mc);

    // The vanilla call is exactly the asset-or-nothing minus K cash-or-nothings.
    const double vanilla = black_scholes_call(S, K, r, sigma, T);
    std::printf("Decomposition  (asset-or-nothing) - K*(cash-or-nothing) = vanilla call\n");
    std::printf("  %.4f - %.0f * %.4f = %.4f   (Black-Scholes %.4f)\n",
                asset_c, K, cash_c, asset_c - K * cash_c, vanilla);

    // Cash-or-nothing parity: a $1 binary call plus a $1 binary put is a sure $1.
    const double cash_p = digital_price(OptionType::Put, DigitalType::CashOrNothing, S, K, r, sigma, T);
    std::printf("\nParity  cash call + cash put = e^-rT:  %.4f + %.4f = %.4f  (e^-rT %.4f)\n",
                cash_c, cash_p, cash_c + cash_p, std::exp(-r * T));
    return 0;
}

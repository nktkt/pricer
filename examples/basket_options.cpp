// basket_options.cpp — multi-asset options on correlated GBM: baskets and spreads.
//
// A basket option pays on a weighted portfolio of assets; a spread option pays on
// the difference of two. Pricing them means simulating correlated GBM (normals
// correlated through the Cholesky factor of the correlation matrix). This demo
// prices a three-asset basket two ways (the exact geometric closed form vs. Monte
// Carlo, and the arithmetic basket by MC), and a two-asset spread three ways
// (Margrabe's exact exchange formula, Kirk's approximation, and MC).
#include "pricer/basket.hpp"
#include "pricer/black_scholes.hpp"

#include <cstdio>
#include <vector>

using namespace pricer;

int main() {
    const double r = 0.05, T = 1.0;
    const long n_paths = 500'000;

    // --- Three-asset basket ---
    const std::vector<double> S{100, 95, 105}, w{1.0 / 3, 1.0 / 3, 1.0 / 3}, sig{0.20, 0.25, 0.30};
    const std::vector<std::vector<double>> corr{
        {1.0, 0.5, 0.3}, {0.5, 1.0, 0.4}, {0.3, 0.4, 1.0}};
    const double K = 100.0;

    std::printf("Basket call on 3 correlated assets   K=%.0f r=%.2f T=%.1f   (%ld paths)\n\n",
                K, r, T, n_paths);
    const double geo_cf = geometric_basket_price(OptionType::Call, S, w, K, r, sig, corr, T);
    const double geo_mc =
        basket_price_mc(OptionType::Call, AverageType::Geometric, S, w, K, r, sig, corr, T, n_paths);
    const double ari_mc =
        basket_price_mc(OptionType::Call, AverageType::Arithmetic, S, w, K, r, sig, corr, T, n_paths);
    std::printf("  geometric basket (closed form): %.4f\n", geo_cf);
    std::printf("  geometric basket (Monte Carlo): %.4f   <- matches closed form\n", geo_mc);
    std::printf("  arithmetic basket (Monte Carlo): %.4f  <- the real product (>= geometric)\n\n",
                ari_mc);

    // --- Two-asset spread / exchange ---
    const double S1 = 100, S2 = 100, sig1 = 0.20, sig2 = 0.30, rho = 0.4;
    std::printf("Spread / exchange call   S1=%.0f S2=%.0f sigma1=%.2f sigma2=%.2f rho=%.1f\n",
                S1, S2, sig1, sig2, rho);
    const double marg = margrabe_exchange_price(S1, S2, sig1, sig2, rho, T);
    const double exch_mc = spread_price_mc(OptionType::Call, S1, S2, 0.0, r, sig1, sig2, rho, T, n_paths);
    std::printf("  exchange (K=0): Margrabe %.4f   Monte Carlo %.4f\n", marg, exch_mc);
    for (double k : {0.0, 5.0, 10.0}) {
        const double kirk = spread_kirk_price(OptionType::Call, S1, S2, k, r, sig1, sig2, rho, T);
        const double mc = spread_price_mc(OptionType::Call, S1, S2, k, r, sig1, sig2, rho, T, n_paths);
        std::printf("  spread K=%-5.1f : Kirk %.4f   Monte Carlo %.4f\n", k, kirk, mc);
    }
    return 0;
}

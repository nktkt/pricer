// vol_surface_demo.cpp — set up market data and use it.
// Walks through three building blocks of a pricing desk: a discount curve, an
// implied-volatility surface back-solved from synthetic market call prices, and
// a least-squares calibration of a single-expiry volatility smile.
#include "pricer/black_scholes.hpp"
#include "pricer/curve.hpp"
#include "pricer/vol_surface.hpp"
#include "pricer/smile.hpp"
#include <cmath>
#include <cstdio>
#include <vector>

using namespace pricer;

int main() {
    // -------------------------------------------------------------------------
    // 1) Discount curve — the term structure we discount cash flows against.
    // -------------------------------------------------------------------------
    const std::vector<double> times = {0.5, 1.0, 2.0, 5.0};
    const std::vector<double> zeros = {0.03, 0.035, 0.04, 0.045};
    DiscountCurve curve(times, zeros);

    std::printf("=== 1) Discount curve ===\n");
    std::printf("%8s | %12s | %12s\n", "t", "zero_rate", "df");
    std::printf("---------|--------------|-------------\n");
    for (double t : {0.5, 1.0, 2.0, 3.0, 5.0})
        std::printf("%8.2f | %12.6f | %12.6f\n", t, curve.zero_rate(t), curve.df(t));
    std::printf("forward_rate(1,2) = %.6f\n\n", curve.forward_rate(1.0, 2.0));

    // -------------------------------------------------------------------------
    // 2) Vol surface — back out implied vols from synthetic market prices.
    // -------------------------------------------------------------------------
    const double S = 100.0, r = 0.05;
    const std::vector<double> expiries = {0.25, 1.0, 2.0};
    const std::vector<double> strikes  = {80, 90, 100, 110, 120};

    // A "true" skew: lower strikes carry higher vol (downside protection bid).
    auto true_vol = [](double K) { return 0.20 + 0.0015 * (100.0 - K); };

    // Quote the market: a Black–Scholes call price per (expiry, strike) node.
    std::vector<std::vector<double>> call_prices(expiries.size(),
                                                 std::vector<double>(strikes.size()));
    for (std::size_t i = 0; i < expiries.size(); ++i)
        for (std::size_t j = 0; j < strikes.size(); ++j)
            call_prices[i][j] =
                black_scholes_call(S, strikes[j], r, true_vol(strikes[j]), expiries[i]);

    // The surface back-solves the implied vol at each node from those prices.
    VolSurface surf(S, r, expiries, strikes, call_prices);

    std::printf("=== 2) Vol surface (recovered implied-vol grid) ===\n");
    std::printf("%8s |", "T \\ K");
    for (double K : strikes) std::printf(" %9.0f", K);
    std::printf("\n---------|");
    for (std::size_t j = 0; j < strikes.size(); ++j) std::printf("----------");
    std::printf("\n");
    for (double T : expiries) {
        std::printf("%8.2f |", T);
        for (double K : strikes) std::printf(" %9.6f", surf.vol(K, T));
        std::printf("\n");
    }

    // Off-node query: interpolated vol and the price it implies.
    std::printf("\nsurf.vol(95, 1.5)        = %.6f\n", surf.vol(95.0, 1.5));
    std::printf("surf.call_price(95, 1.5) = %.6f\n", surf.call_price(95.0, 1.5));

    // A node should reproduce the true vol it was built from.
    const double node_recovered = surf.vol(80.0, 1.0);
    const double node_true = true_vol(80.0);
    std::printf("node check: surf.vol(80,1.0) = %.6f vs true = %.6f (diff %.2e)\n\n",
                node_recovered, node_true, std::fabs(node_recovered - node_true));

    // -------------------------------------------------------------------------
    // 3) Smile calibration — fit a quadratic in log-moneyness to the T=1 slice.
    // -------------------------------------------------------------------------
    std::vector<double> vols_1y(strikes.size());
    for (std::size_t j = 0; j < strikes.size(); ++j)
        vols_1y[j] = surf.vol(strikes[j], 1.0);

    const Smile smile = calibrate_smile(S, strikes, vols_1y);
    std::printf("=== 3) Smile calibration (T = 1.0) ===\n");
    std::printf("fitted vol(k) = a + b*k + c*k^2,  k = log(K/S)\n");
    std::printf("a = %.6f, b = %.6f, c = %.6f\n", smile.a, smile.b, smile.c);

    double max_err = 0.0;
    for (std::size_t j = 0; j < strikes.size(); ++j)
        max_err = std::fmax(max_err, std::fabs(smile.vol(strikes[j]) - vols_1y[j]));
    std::printf("max abs reproduction error over strikes = %.2e\n", max_err);

    return 0;
}

// xva_demo.cpp — exposure simulation and CVA/DVA for a single trade.
//
// Prices the counterparty credit risk on a one-year European call. The steps are
// the standard xVA pipeline:
//   1. Simulate the underlying forward under risk-neutral GBM (the scenario
//      engine) and revalue the option at each monthly grid point, giving the
//      expected positive/negative exposure profile.
//   2. Turn a flat credit spread into a survival (default-probability) curve.
//   3. Integrate exposure × marginal default probability × discount factor into
//      CVA (counterparty default) and DVA (own default).
#include "pricer/black_scholes.hpp"
#include "pricer/curve.hpp"
#include "pricer/xva.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <vector>

using namespace pricer;

int main() {
    const double S = 100, K = 100, r = 0.03, sigma = 0.20, T = 1.0;
    const double cp_spread = 0.020, own_spread = 0.010, recovery = 0.40;
    const long n_paths = 400'000;

    // Monthly exposure grid out to expiry.
    std::vector<double> grid;
    for (int i = 1; i <= 12; ++i) grid.push_back(i / 12.0);

    // Mark-to-market of a long European call at time t given spot S_t.
    auto call_mtm = [&](double t, double St) {
        return (T - t > 1e-8) ? black_scholes_call(St, K, r, sigma, T - t) : std::max(St - K, 0.0);
    };

    const ExposureProfile ep = exposure_profile_gbm(call_mtm, S, r, sigma, grid, n_paths);

    const DiscountCurve disc({0.5, 1.0, 2.0}, {r, r, r});  // flat at r
    const SurvivalCurve cp = SurvivalCurve::from_spread(cp_spread, recovery);
    const SurvivalCurve own = SurvivalCurve::from_spread(own_spread, recovery);

    const double price = black_scholes_call(S, K, r, sigma, T);
    std::printf("European call  S=%.0f K=%.0f r=%.2f sigma=%.2f T=%.1f\n", S, K, r, sigma, T);
    std::printf("Risk-free price = %.4f   (%ld paths, monthly grid)\n\n", price, n_paths);

    std::printf("Expected exposure profile (EPE; ENE is ~0 for a long option):\n");
    std::printf("  %-8s  %12s  %12s\n", "t(yr)", "EPE", "disc.EPE");
    for (std::size_t i = 0; i < ep.times.size(); ++i)
        std::printf("  %-8.4f  %12.5f  %12.5f\n", ep.times[i], ep.epe[i],
                    disc.df(ep.times[i]) * ep.epe[i]);

    const double cva_v = cva(ep, cp, disc, recovery);
    const double dva_v = dva(ep, own, disc, recovery);
    const double bcva_v = bcva(ep, cp, own, disc, recovery, recovery);

    std::printf("\nCredit assumptions: counterparty spread=%.0fbp, own spread=%.0fbp, recovery=%.0f%%\n",
                cp_spread * 1e4, own_spread * 1e4, recovery * 100);
    std::printf("CVA  = %.5f   (%.2f%% of price)\n", cva_v, 100 * cva_v / price);
    std::printf("DVA  = %.5f\n", dva_v);
    std::printf("BCVA = CVA - DVA = %.5f\n", bcva_v);
    std::printf("Credit-adjusted value to us = price - BCVA = %.5f\n", price - bcva_v);
    return 0;
}

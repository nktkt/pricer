// Tests for counterparty valuation adjustments (xva.hpp).
#include "check.hpp"
#include "pricer/black_scholes.hpp"
#include "pricer/curve.hpp"
#include "pricer/xva.hpp"

#include <algorithm>
#include <cmath>
#include <vector>

using namespace pricer;

int main() {
    // --- survival curve basics ---
    SurvivalCurve sc{0.03};
    check::approx("survival(0) = 1", sc.survival(0.0), 1.0, 1e-12);
    check::is_true("survival decreasing", sc.survival(2.0) < sc.survival(1.0));
    check::is_true("default prob >= 0", sc.default_prob(1.0, 2.0) > 0.0);
    const SurvivalCurve from_spread = SurvivalCurve::from_spread(0.012, 0.4);
    check::approx("hazard from spread", from_spread.hazard, 0.012 / 0.6, 1e-12);

    const double S = 100, K = 100, r = 0.03, sigma = 0.20, T = 1.0;
    const DiscountCurve disc({0.25, 0.5, 1.0, 2.0}, {r, r, r, r});  // flat at r
    std::vector<double> grid;
    for (int i = 1; i <= 12; ++i) grid.push_back(i / 12.0);  // monthly to expiry

    // --- a long European call: MtM is never negative ---
    auto call_mtm = [&](double t, double St) {
        return (T - t > 1e-8) ? black_scholes_call(St, K, r, sigma, T - t) : std::max(St - K, 0.0);
    };
    const ExposureProfile ep = exposure_profile_gbm(call_mtm, S, r, sigma, grid, 200'000);
    double max_ene = 0.0;
    for (double e : ep.ene) max_ene = std::max(max_ene, e);
    check::approx("long call has no negative exposure", max_ene, 0.0, 1e-12);

    const SurvivalCurve no_default{0.0};
    check::approx("CVA = 0 when hazard = 0", cva(ep, no_default, disc, 0.4), 0.0, 1e-12);

    const SurvivalCurve cp = SurvivalCurve::from_spread(0.02, 0.4);
    const double cva_base = cva(ep, cp, disc, 0.4);
    check::is_true("CVA > 0 with credit risk", cva_base > 0.0);
    check::approx("DVA = 0 for a long call", dva(ep, cp, disc, 0.4), 0.0, 1e-12);
    check::approx("CVA = 0 at full recovery", cva(ep, cp, disc, 1.0), 0.0, 1e-12);

    const SurvivalCurve wider = SurvivalCurve::from_spread(0.04, 0.4);
    check::is_true("CVA increases with spread", cva(ep, wider, disc, 0.4) > cva_base);

    // --- a forward contract: MtM swings both ways, so DVA is non-zero too ---
    auto fwd_mtm = [&](double t, double St) { return St - K * std::exp(-r * (T - t)); };
    const ExposureProfile epf = exposure_profile_gbm(fwd_mtm, S, r, sigma, grid, 200'000);
    const double cva_f = cva(epf, cp, disc, 0.4), dva_f = dva(epf, cp, disc, 0.4);
    check::is_true("forward CVA > 0", cva_f > 0.0);
    check::is_true("forward DVA > 0", dva_f > 0.0);
    check::approx("BCVA = CVA - DVA", bcva(epf, cp, cp, disc, 0.4, 0.4), cva_f - dva_f, 1e-12);

    return check::report("xva");
}

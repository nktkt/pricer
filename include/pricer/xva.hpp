// pricer/xva.hpp — counterparty valuation adjustments (CVA / DVA) and the
// exposure-simulation scenario engine they rest on.
//
// A derivative's value to you depends on the counterparty actually paying. CVA
// is the market price of that credit risk: the expected discounted loss if the
// counterparty defaults while the trade is in your favour. DVA is the mirror
// image for your own default. Both are driven by the *exposure profile* — how
// the position's mark-to-market is expected to evolve over its life — which we
// obtain by simulating market scenarios forward and revaluing at each step.
//
//   CVA = (1 - R_cp ) * Σ DF(t_i) * EPE(t_i) * PD_cp (t_{i-1}, t_i)
//   DVA = (1 - R_own) * Σ DF(t_i) * ENE(t_i) * PD_own(t_{i-1}, t_i)
//
// where EPE/ENE are the expected positive/negative exposures, DF the discount
// factor, R the recovery rate and PD the marginal default probability.
#pragma once
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <vector>

#include "pricer/curve.hpp"  // DiscountCurve
#include "pricer/rng.hpp"    // cb_normal

namespace pricer {

// Constant-hazard survival model: S(t) = exp(-h t). The marginal default
// probability over (a, b] is S(a) - S(b). A flat credit spread implies a hazard
// of spread / (1 - recovery).
struct SurvivalCurve {
    double hazard = 0.0;

    double survival(double t) const { return std::exp(-hazard * t); }
    double default_prob(double a, double b) const { return survival(a) - survival(b); }

    static SurvivalCurve from_spread(double spread, double recovery) {
        const double lgd = 1.0 - recovery;
        return {spread / (lgd > 1e-12 ? lgd : 1e-12)};
    }
};

// Expected exposure profile on a time grid (undiscounted; discounting is applied
// in cva()/dva()). EPE = E[max(V,0)], ENE = E[max(-V,0)] at each grid time.
struct ExposureProfile {
    std::vector<double> times;
    std::vector<double> epe;
    std::vector<double> ene;
};

// Simulate the exposure profile of a position under risk-neutral GBM. `mtm(t, S)`
// returns the position's mark-to-market value at time t given spot S. Paths use
// the counter-based RNG, so the profile is reproducible for a given seed.
template <class Mtm>
ExposureProfile exposure_profile_gbm(Mtm mtm, double S0, double r, double sigma,
                                     const std::vector<double>& grid, long n_paths,
                                     std::uint64_t seed = 12345) {
    if (grid.empty()) throw std::invalid_argument("exposure_profile_gbm: empty time grid");
    if (n_paths <= 0) throw std::invalid_argument("exposure_profile_gbm: n_paths must be > 0");

    const std::size_t m = grid.size();
    std::vector<double> sum_pos(m, 0.0), sum_neg(m, 0.0);

    for (long p = 0; p < n_paths; ++p) {
        double S = S0, t_prev = 0.0;
        for (std::size_t i = 0; i < m; ++i) {
            const double dt = grid[i] - t_prev;
            if (dt < 0.0) throw std::invalid_argument("exposure_profile_gbm: grid not increasing");
            // One independent draw per (path, step); distinct counters -> independent.
            const std::uint64_t counter = static_cast<std::uint64_t>(p) * m + i;
            const double z = cb_normal(seed, counter);
            S *= std::exp((r - 0.5 * sigma * sigma) * dt + sigma * std::sqrt(dt) * z);
            const double V = mtm(grid[i], S);
            sum_pos[i] += V > 0.0 ? V : 0.0;
            sum_neg[i] += V < 0.0 ? -V : 0.0;
            t_prev = grid[i];
        }
    }

    ExposureProfile ep;
    ep.times = grid;
    ep.epe.resize(m);
    ep.ene.resize(m);
    const double inv_n = 1.0 / static_cast<double>(n_paths);
    for (std::size_t i = 0; i < m; ++i) {
        ep.epe[i] = sum_pos[i] * inv_n;
        ep.ene[i] = sum_neg[i] * inv_n;
    }
    return ep;
}

// Credit valuation adjustment: expected discounted loss from the counterparty's
// default over the life of the trade.
inline double cva(const ExposureProfile& ep, const SurvivalCurve& counterparty,
                  const DiscountCurve& disc, double recovery) {
    double acc = 0.0, t_prev = 0.0;
    for (std::size_t i = 0; i < ep.times.size(); ++i) {
        const double pd = counterparty.default_prob(t_prev, ep.times[i]);
        acc += disc.df(ep.times[i]) * ep.epe[i] * pd;
        t_prev = ep.times[i];
    }
    return (1.0 - recovery) * acc;
}

// Debit valuation adjustment: the mirror of CVA for your own default (uses the
// expected negative exposure and your own survival curve).
inline double dva(const ExposureProfile& ep, const SurvivalCurve& own, const DiscountCurve& disc,
                  double recovery_own) {
    double acc = 0.0, t_prev = 0.0;
    for (std::size_t i = 0; i < ep.times.size(); ++i) {
        const double pd = own.default_prob(t_prev, ep.times[i]);
        acc += disc.df(ep.times[i]) * ep.ene[i] * pd;
        t_prev = ep.times[i];
    }
    return (1.0 - recovery_own) * acc;
}

// Bilateral CVA: the net credit adjustment to the risk-free value, BCVA = CVA - DVA.
inline double bcva(const ExposureProfile& ep, const SurvivalCurve& counterparty,
                   const SurvivalCurve& own, const DiscountCurve& disc, double recovery_cp,
                   double recovery_own) {
    return cva(ep, counterparty, disc, recovery_cp) - dva(ep, own, disc, recovery_own);
}

}  // namespace pricer

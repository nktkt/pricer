// mc_greeks.cpp — Monte Carlo Greeks vs. closed-form Black–Scholes Greeks.
//
// Shows two ways to get sensitivities from a simulation:
//   * bump-and-revalue with common random numbers (reprice at S±h / sigma±h
//     using the same seed, so the noise cancels), and
//   * a pathwise derivative for the call's delta (the adjoint-AD-style estimator,
//     no bumping).
// Both are checked against the analytic Greeks.
#include "pricer/black_scholes.hpp"
#include "pricer/greeks_mc.hpp"
#include <cstdio>

using namespace pricer;

int main() {
    const double S = 100, K = 100, r = 0.05, sigma = 0.20, T = 1.0;
    const long n = 8'000'000;
    auto call = [K](double ST) { return ST > K ? ST - K : 0.0; };

    const Greeks g = black_scholes_greeks(OptionType::Call, S, K, r, sigma, T);
    const double d_crn = mc::mc_delta_crn(call, S, r, sigma, T, n);
    const double v_crn = mc::mc_vega_crn(call, S, r, sigma, T, n);
    const double d_pw  = mc::mc_call_delta_pathwise(S, K, r, sigma, T, n);

    std::printf("call: S=%.0f K=%.0f r=%.2f sigma=%.2f T=%.0f  paths=%ld\n\n", S, K, r, sigma, T, n);
    std::printf("%-22s | %12s | %12s\n", "greek", "Monte Carlo", "closed form");
    std::printf("-----------------------|--------------|--------------\n");
    std::printf("%-22s | %12.6f | %12.6f\n", "delta (bump + CRN)", d_crn, g.delta);
    std::printf("%-22s | %12.6f | %12.6f\n", "delta (pathwise)",   d_pw,  g.delta);
    std::printf("%-22s | %12.6f | %12.6f\n", "vega  (bump + CRN)", v_crn, g.vega);
    std::printf("\nCommon random numbers make the bump estimates stable; the pathwise\n"
                "delta needs no bump at all (this is what adjoint AD computes).\n");
    return 0;
}

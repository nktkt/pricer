// Tests for Monte Carlo Greeks against closed-form Black–Scholes Greeks.
#include "check.hpp"
#include "pricer/black_scholes.hpp"
#include "pricer/greeks_mc.hpp"

using namespace pricer;

int main() {
    const double S = 100, K = 100, r = 0.05, sigma = 0.20, T = 1.0;
    auto call = [K](double ST) { return ST > K ? ST - K : 0.0; };
    const Greeks g = black_scholes_greeks(OptionType::Call, S, K, r, sigma, T);
    const long n = 4'000'000;

    // Common-random-number bump greeks have low variance, so a tight tolerance.
    check::approx("delta (bump+CRN)", mc::mc_delta_crn(call, S, r, sigma, T, n), g.delta, 5e-3);
    check::approx("vega (bump+CRN)",  mc::mc_vega_crn(call, S, r, sigma, T, n),  g.vega,  5e-2);

    // Pathwise delta (the adjoint-AD-style estimator).
    check::approx("delta (pathwise)", mc::mc_call_delta_pathwise(S, K, r, sigma, T, n),
                  g.delta, 5e-3);

    return check::report("greeks_mc");
}

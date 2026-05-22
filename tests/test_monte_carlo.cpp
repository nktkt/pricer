// Tests that the generic Monte Carlo engine converges to the analytic price.
#include "check.hpp"
#include "pricer/black_scholes.hpp"
#include "pricer/monte_carlo.hpp"

using namespace pricer;

int main() {
    const double S = 100, K = 100, r = 0.05, sigma = 0.20, T = 1.0;
    const long n = 4'000'000;

    // Call payoff supplied as a lambda — the engine stays instrument-agnostic.
    const double mc_call = mc::price_terminal(
        [K](double ST) { return ST > K ? ST - K : 0.0; }, S, r, sigma, T, n, /*seed=*/777);
    const double mc_put = mc::price_terminal(
        [K](double ST) { return ST < K ? K - ST : 0.0; }, S, r, sigma, T, n, /*seed=*/778);

    check::approx("MC call vs BS", mc_call, black_scholes_call(S, K, r, sigma, T), 0.02);
    check::approx("MC put vs BS",  mc_put,  black_scholes_put(S, K, r, sigma, T),  0.02);

    // Determinism: same seed must give the exact same result.
    const double a = mc::price_terminal([K](double ST){ return ST > K ? ST - K : 0.0; }, S, r, sigma, T, 100000, 42);
    const double b = mc::price_terminal([K](double ST){ return ST > K ? ST - K : 0.0; }, S, r, sigma, T, 100000, 42);
    check::is_true("seeded result is reproducible", a == b);

    return check::report("monte_carlo");
}

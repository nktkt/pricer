// Tests for the counter-based RNG and its terminal-value Monte Carlo.
#include "check.hpp"
#include "pricer/black_scholes.hpp"
#include "pricer/rng.hpp"
#include <cmath>

using namespace pricer;

int main() {
    // Stateless & reproducible: the same (seed, counter) always gives the same draw.
    check::is_true("uniform reproducible", cb_uniform(7, 42) == cb_uniform(7, 42));
    check::is_true("normal reproducible", cb_normal(7, 42) == cb_normal(7, 42));
    check::is_true("different counters differ", cb_uniform(7, 1) != cb_uniform(7, 2));
    check::is_true("different seeds differ", cb_uniform(1, 5) != cb_uniform(2, 5));

    // Uniforms stay strictly inside (0,1); rough mean ~ 0.5 over a large sample.
    double su = 0.0;
    const long N = 2'000'000;
    bool in_range = true;
    for (long i = 0; i < N; ++i) {
        const double u = cb_uniform(99, static_cast<std::uint64_t>(i));
        if (u <= 0.0 || u >= 1.0) in_range = false;
        su += u;
    }
    check::is_true("uniforms in (0,1)", in_range);
    check::approx("uniform mean ~ 0.5", su / N, 0.5, 2e-3);

    // Counter-based MC converges to Black–Scholes.
    const double S = 100, K = 100, r = 0.05, sigma = 0.20, T = 1.0;
    auto call = [K](double ST) { return ST > K ? ST - K : 0.0; };
    const double mc = mc::price_terminal_cb(call, S, r, sigma, T, 4'000'000);
    check::approx("counter-based MC vs BS", mc, black_scholes_call(S, K, r, sigma, T), 0.02);

    // The whole point: re-running gives the identical price (no hidden state).
    const double mc2 = mc::price_terminal_cb(call, S, r, sigma, T, 4'000'000);
    check::is_true("MC reproducible run-to-run", mc == mc2);

    return check::report("rng");
}

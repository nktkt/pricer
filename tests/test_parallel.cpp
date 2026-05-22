// Tests for the deterministic parallel Monte Carlo engine.
#include "check.hpp"
#include "pricer/black_scholes.hpp"
#include "pricer/parallel.hpp"

using namespace pricer;

int main() {
    const double S = 100, K = 100, r = 0.05, sigma = 0.20, T = 1.0;
    auto call = [K](double ST) { return ST > K ? ST - K : 0.0; };
    const double bs = black_scholes_call(S, K, r, sigma, T);

    const long n = 4'000'000;
    const double p1 = mc::price_terminal_parallel(call, S, r, sigma, T, n, 123, 1);
    const double p4 = mc::price_terminal_parallel(call, S, r, sigma, T, n, 123, 4);
    const double p8 = mc::price_terminal_parallel(call, S, r, sigma, T, n, 123, 8);

    // The whole point: identical result no matter how many threads ran it.
    check::is_true("deterministic: 1 vs 4 threads", p1 == p4);
    check::is_true("deterministic: 1 vs 8 threads", p1 == p8);
    check::approx("parallel price vs BS", p1, bs, 0.02);

    return check::report("parallel");
}

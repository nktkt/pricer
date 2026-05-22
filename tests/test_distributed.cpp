// Tests that sharded Monte Carlo is reproducible across any shard count.
#include "check.hpp"
#include "pricer/black_scholes.hpp"
#include "pricer/distributed.hpp"

using namespace pricer;

int main() {
    const double S = 100, K = 100, r = 0.05, sigma = 0.20, T = 1.0;
    auto call = [K](double ST) { return ST > K ? ST - K : 0.0; };
    const long n = 4'000'000;

    const double p1 = mc::price_terminal_sharded(call, S, r, sigma, T, n, 1);
    const double p3 = mc::price_terminal_sharded(call, S, r, sigma, T, n, 3);
    const double p8 = mc::price_terminal_sharded(call, S, r, sigma, T, n, 8);

    // The whole point: bit-for-bit identical regardless of how blocks are sharded.
    check::is_true("bitwise identical: 1 vs 3 shards", p1 == p3);
    check::is_true("bitwise identical: 1 vs 8 shards", p1 == p8);
    check::approx("sharded price vs BS", p1, black_scholes_call(S, K, r, sigma, T), 0.02);

    return check::report("distributed");
}

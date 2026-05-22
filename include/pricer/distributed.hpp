// pricer/distributed.hpp — building blocks for distributed Monte Carlo.
//
// The work is a FIXED B-block decomposition: block b always uses the same
// sub-seed and the same share of paths, and its payoff sum is computed
// independently. Whoever computes a block gets the same number, and the final
// reduction always sums the B block-sums in canonical order 0..B-1. Therefore
// the result is BIT-FOR-BIT identical no matter how the blocks are distributed
// across shards / threads / processes / nodes — the property a distributed MC
// needs for reproducible risk. (A multi-process driver using these blocks lives
// in examples/distributed_mc.cpp.)
#pragma once
#include <cmath>
#include <cstdint>
#include <random>
#include <vector>

#include "pricer/parallel.hpp"  // detail::splitmix64

namespace pricer::mc {

// Compute the payoff sum of each block in [b0, b1) of a B-block decomposition of
// `n_paths`, writing them to out[0 .. (b1-b0)). Deterministic per block.
template <class Payoff>
void fill_block_sums(Payoff payoff, double S, double r, double sigma, double T, long n_paths,
                     int B, int b0, int b1, std::uint64_t seed, double* out) {
    const double drift = (r - 0.5 * sigma * sigma) * T;
    const double vol = sigma * std::sqrt(T);
    const long base = n_paths / B, rem = n_paths % B;
    for (int b = b0; b < b1; ++b) {
        const long cnt = base + (b < rem ? 1 : 0);
        std::mt19937_64 rng(detail::splitmix64(seed + static_cast<std::uint64_t>(b) + 1));
        std::normal_distribution<double> Z(0.0, 1.0);
        double s = 0.0;
        for (long i = 0; i < cnt; ++i) s += payoff(S * std::exp(drift + vol * Z(rng)));
        out[b - b0] = s;
    }
}

// Discount an aggregated payoff sum to a price.
inline double aggregate_price(double total_sum, long n_paths, double r, double T) {
    return std::exp(-r * T) * (total_sum / static_cast<double>(n_paths));
}

// In-process reference: split B blocks across `shards` contiguous groups, then
// reduce all B block-sums in canonical order. The returned price is identical
// for any `shards` — the same invariant the multi-process driver relies on.
template <class Payoff>
double price_terminal_sharded(Payoff payoff, double S, double r, double sigma, double T,
                              long n_paths, int shards, int B = 512, std::uint64_t seed = 12345) {
    std::vector<double> block(B, 0.0);
    for (int s = 0; s < shards; ++s) {
        const int b0 = static_cast<int>(static_cast<long>(B) * s / shards);
        const int b1 = static_cast<int>(static_cast<long>(B) * (s + 1) / shards);
        fill_block_sums(payoff, S, r, sigma, T, n_paths, B, b0, b1, seed, &block[b0]);
    }
    double total = 0.0;
    for (int b = 0; b < B; ++b) total += block[b];  // canonical-order reduction
    return aggregate_price(total, n_paths, r, T);
}

}  // namespace pricer::mc

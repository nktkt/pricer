// pricer/simd_mc.hpp — SIMD-vectorized Monte Carlo path generation.
//
// The counter-based RNG (rng.hpp) makes draw i a pure function of (seed, i), so a
// batch of W counters can be processed in one SIMD register with no carried
// state. This header generates W GBM terminal prices at a time — vectorized
// integer mixing, uniform conversion, inverse-normal CDF and exp — then applies
// the (scalar) payoff per lane. The expensive path-generation work runs on the
// SIMD units; the result agrees with the scalar engine to floating-point
// tolerance (the integer mixing is bit-identical; the transcendentals match to
// ~1e-11). Falls back to the scalar engine where the vector extensions are
// unavailable.
#pragma once
#include <cmath>
#include <cstdint>

#include "pricer/rng.hpp"
#include "pricer/simd.hpp"

namespace pricer::mc {

#if PRICER_HAVE_SIMD

// W uniforms in (0,1) for counters [base, base+W). Bit-identical to W scalar
// cb_uniform calls: same SplitMix64 mixing, same 53-bit mantissa scaling.
inline simd::vd cb_uniform_simd(std::uint64_t seed, std::uint64_t base) {
    simd::vu counter = simd::splat_u(base) + simd::vu{0, 1, 2, 3};
    simd::vu x = simd::splat_u(seed) + counter * simd::splat_u(0x9E3779B97F4A7C15ull);
    // SplitMix64 (matches mc::detail::splitmix64), vectorized.
    x = x + simd::splat_u(0x9E3779B97F4A7C15ull);
    x = (x ^ (x >> 30)) * simd::splat_u(0xBF58476D1CE4E5B9ull);
    x = (x ^ (x >> 27)) * simd::splat_u(0x94D049BB133111EBull);
    x = x ^ (x >> 31);
    simd::vd z = __builtin_convertvector(x >> 11, simd::vd);
    return (z + simd::splat(0.5)) * simd::splat(1.0 / 9007199254740992.0);
}

// W standard-normal draws for counters [base, base+W).
inline simd::vd cb_normal_simd(std::uint64_t seed, std::uint64_t base) {
    return simd::inv_norm_cdf(cb_uniform_simd(seed, base));
}

// Terminal-value Monte Carlo with SIMD path generation. Same contract and result
// (to tolerance) as price_terminal_cb, but W paths are generated per iteration.
// The leftover paths (n mod W) use the scalar kernel.
template <class Payoff>
double price_terminal_cb_simd(Payoff payoff, double S, double r, double sigma, double T, long n,
                              std::uint64_t seed = 12345) {
    const double drift = (r - 0.5 * sigma * sigma) * T;
    const double vol = sigma * std::sqrt(T);
    const simd::vd vS = simd::splat(S), vdrift = simd::splat(drift), vvol = simd::splat(vol);

    double sum = 0.0;
    long i = 0;
    const long w = simd::kWidth;
    for (; i + w <= n; i += w) {
        simd::vd z = cb_normal_simd(seed, static_cast<std::uint64_t>(i));
        simd::vd ST = vS * simd::v_exp(vdrift + vvol * z);
        for (int j = 0; j < w; ++j) sum += payoff(ST[j]);  // scalar payoff per lane
    }
    for (; i < n; ++i)  // scalar tail for the final n mod W paths
        sum += payoff(S * std::exp(drift + vol * cb_normal(seed, static_cast<std::uint64_t>(i))));

    return std::exp(-r * T) * (sum / static_cast<double>(n));
}

#else  // no vector extensions: fall back to the scalar counter-based engine.

template <class Payoff>
double price_terminal_cb_simd(Payoff payoff, double S, double r, double sigma, double T, long n,
                              std::uint64_t seed = 12345) {
    return price_terminal_cb(payoff, S, r, sigma, T, n, seed);
}

#endif  // PRICER_HAVE_SIMD

}  // namespace pricer::mc

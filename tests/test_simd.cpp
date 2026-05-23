// Tests for the SIMD layer (simd.hpp) and SIMD path generation (simd_mc.hpp).
#include "check.hpp"
#include "pricer/black_scholes.hpp"
#include "pricer/qmc.hpp"  // scalar inv_norm_cdf
#include "pricer/rng.hpp"  // scalar cb_uniform / cb_normal / price_terminal_cb
#include "pricer/simd.hpp"
#include "pricer/simd_mc.hpp"

#include <cmath>
#include <cstdint>

using namespace pricer;

int main() {
#if PRICER_HAVE_SIMD
    // --- vectorized transcendentals vs the standard library ---
    double max_exp = 0, max_log = 0, max_sqrt = 0;
    for (double x = -30; x <= 30; x += 0.05) {
        const double e = std::fabs(simd::v_exp(simd::splat(x))[0] - std::exp(x)) / std::exp(x);
        if (e > max_exp) max_exp = e;
    }
    for (double x = 1e-4; x <= 1e4; x *= 1.01) {
        const double e = std::fabs(simd::v_log(simd::splat(x))[0] - std::log(x)) /
                         (std::fabs(std::log(x)) + 1e-12);
        if (e > max_log) max_log = e;
    }
    for (double x = 0; x <= 1e4; x += 0.7) {
        const double e =
            std::fabs(simd::v_sqrt(simd::splat(x))[0] - std::sqrt(x)) / (std::sqrt(x) + 1e-12);
        if (e > max_sqrt) max_sqrt = e;
    }
    check::is_true("v_exp  accurate (<1e-9)", max_exp < 1e-9);
    check::is_true("v_log  accurate (<1e-9)", max_log < 1e-9);
    check::is_true("v_sqrt accurate (<1e-9)", max_sqrt < 1e-9);

    // --- vectorized inverse-normal CDF vs the scalar version ---
    double max_inv = 0;
    for (double p = 1e-5; p < 1.0; p += 1e-4) {
        const double e = std::fabs(simd::inv_norm_cdf(simd::splat(p))[0] - inv_norm_cdf(p));
        if (e > max_inv) max_inv = e;
    }
    check::is_true("simd inv_norm_cdf vs scalar (<1e-7)", max_inv < 1e-7);

    // --- SIMD uniforms are bit-identical to the scalar counter-based RNG ---
    bool unif_exact = true;
    for (std::uint64_t base = 0; base < 40; base += simd::kWidth) {
        const simd::vd u = mc::cb_uniform_simd(123, base);
        for (int j = 0; j < simd::kWidth; ++j)
            if (u[j] != cb_uniform(123, base + static_cast<std::uint64_t>(j))) unif_exact = false;
    }
    check::is_true("simd uniforms bit-identical to scalar", unif_exact);

    // --- SIMD normals match the scalar normals to tolerance ---
    double max_z = 0;
    for (std::uint64_t base = 0; base < 4000; base += simd::kWidth) {
        const simd::vd z = mc::cb_normal_simd(123, base);
        for (int j = 0; j < simd::kWidth; ++j) {
            const double e = std::fabs(z[j] - cb_normal(123, base + static_cast<std::uint64_t>(j)));
            if (e > max_z) max_z = e;
        }
    }
    check::is_true("simd normals match scalar (<1e-7)", max_z < 1e-7);

    // --- SIMD terminal-value MC matches Black–Scholes and the scalar engine ---
    const double S = 100, K = 100, r = 0.05, sigma = 0.20, T = 1.0;
    auto call = [K](double ST) { return ST > K ? ST - K : 0.0; };
    const long n = 4'000'000;
    const double bs = black_scholes_call(S, K, r, sigma, T);
    const double simd_price = mc::price_terminal_cb_simd(call, S, r, sigma, T, n);
    const double scalar_price = mc::price_terminal_cb(call, S, r, sigma, T, n);
    check::approx("simd MC vs Black-Scholes", simd_price, bs, 0.02);
    check::approx("simd MC vs scalar cb MC", simd_price, scalar_price, 1e-3);

    // Same seed/count -> identical price run to run (no hidden state).
    const double simd_price2 = mc::price_terminal_cb_simd(call, S, r, sigma, T, n);
    check::is_true("simd MC reproducible run-to-run", simd_price == simd_price2);

    // A path count not divisible by the lane width still prices correctly.
    const double odd = mc::price_terminal_cb_simd(call, S, r, sigma, T, n + 3);
    check::approx("simd MC handles n mod W tail", odd, bs, 0.02);
#else
    check::is_true("SIMD unavailable: scalar fallback compiles", true);
#endif
    return check::report("simd");
}
